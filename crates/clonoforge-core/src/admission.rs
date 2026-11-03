//! Conservative phase planning for prototype bounded workflow.
//! This is an explicit model, not an exact RSS predictor. It does not replace
//! allocator/OS overhead measurement or enforce parser limits before parsing.
use crate::model::invalid;
use std::io::Result;

const MIB: u64 = 1024 * 1024;
const PARSER_LINE: u64 = 16 * MIB;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Request {
    pub budget: u64,
    pub workers: u64,
    pub window: u64,
    pub merge_fan_in: u64,
    pub max_record_bytes: u64,
    pub max_group_records: u64,
    pub sort_chunk_bytes: u64,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Plan {
    pub fixed: u64,
    pub parser_and_decompressor: u64,
    pub grouping: u64,
    pub worker_pipeline: u64,
    pub output_staging: u64,
    pub safety_margin: u64,
    pub total: u64,
    pub limiting_phase: &'static str,
    pub record_weight: u64,
    pub group_weight: u64,
    pub max_line_bytes: u64,
}
fn mul(a: u64, b: u64) -> Result<u64> {
    a.checked_mul(b)
        .ok_or_else(|| invalid("admission: arithmetic overflow"))
}
fn add(a: u64, b: u64) -> Result<u64> {
    a.checked_add(b)
        .ok_or_else(|| invalid("admission: arithmetic overflow"))
}
impl Request {
    pub fn plan(&self) -> Result<Plan> {
        self.plan_bounded(50)
    }
    #[cfg(feature = "experimental-large")]
    pub fn plan_experimental(&self) -> Result<Plan> {
        if ![100, 250, 500].contains(&self.max_group_records) {
            return Err(invalid("admission: unsupported experimental group limit"));
        }
        self.plan_bounded(500)
    }
    fn plan_bounded(&self, ceiling: u64) -> Result<Plan> {
        if self.workers == 0
            || self.window == 0
            || self.merge_fan_in < 2
            || self.max_record_bytes == 0
            || self.max_record_bytes > PARSER_LINE
            || self.max_group_records == 0
            || self.max_group_records > ceiling
            || self.sort_chunk_bytes == 0
        {
            return Err(invalid("admission: invalid configuration"));
        }
        // Extremely annotation-dense headers can cause many
        // heap-allocated strings and doubling of Vec capacity. 64 bytes per
        // input byte is a deliberately conservative budgeting assumption, not
        // a mathematically portable allocator bound.
        // The dedicated parser checks this record limit before header expansion.
        // A line may include a prefix and CRLF. The legacy 16 MiB cap still applies.
        let max_line_bytes = add(self.max_record_bytes, 3)?.min(PARSER_LINE);
        let parser_and_decompressor = add(mul(max_line_bytes, 64)?, 136 * MIB)?;
        let fixed = 64 * MIB;
        // Same accounting as the existing external sort's per-item estimate:
        // 40× header allowance + 4× sequence + 4096 bytes of metadata.
        // max_record bounds the sum of header and sequence/qualities.
        let record_weight = add(mul(add(self.max_record_bytes, 16)?, 42)?, 4096)?;
        if self.sort_chunk_bytes < record_weight {
            return Err(invalid(
                "admission: sort chunk cannot hold one admitted record",
            ));
        }
        let group_weight = mul(record_weight, self.max_group_records)?;
        // Three live chunks in sort/merge include input, validation and writer;
        // merge heads, file buffers and readback are additionally budgeted.
        let grouping = add(
            mul(self.sort_chunk_bytes, 3)?,
            mul(self.merge_fan_in, add(record_weight, 8192)?)?,
        )?;
        // Each in-flight group may coexist with worker/output serialization
        // copies; worker stacks are explicitly counted at 8 MiB each.
        let worker_pipeline = add(
            mul(mul(self.window, group_weight)?, 3)?,
            mul(self.workers, 8 * MIB)?,
        )?;
        let output_staging = add(mul(group_weight, 2)?, 16 * MIB)?;
        let subtotal = add(
            add(
                add(fixed, parser_and_decompressor)?,
                add(grouping, worker_pipeline)?,
            )?,
            output_staging,
        )?;
        let safety_margin = add(subtotal, 3)? / 4;
        let total = add(subtotal, safety_margin)?;
        if self.budget < total {
            return Err(invalid(&format!(
                "admission: budget {} below modeled minimum {}",
                self.budget, total
            )));
        }
        let limiting_phase = [
            ("parser", parser_and_decompressor),
            ("grouping", grouping),
            ("worker_pipeline", worker_pipeline),
            ("output_staging", output_staging),
        ]
        .into_iter()
        .max_by_key(|x| x.1)
        .unwrap()
        .0;
        Ok(Plan {
            fixed,
            parser_and_decompressor,
            grouping,
            worker_pipeline,
            output_staging,
            safety_margin,
            total,
            limiting_phase,
            record_weight,
            group_weight,
            max_line_bytes,
        })
    }
    /// Choose the largest worker count that fits, preferring larger runs when
    /// admitted and falling back to a smaller run at the same worker count.
    /// This deterministic model is not a hard RSS guarantee.
    pub fn select(
        budget: u64,
        requested_workers: u64,
        max_record_bytes: u64,
        max_group_records: u64,
    ) -> Result<(Self, Plan)> {
        Self::select_bounded(
            budget,
            requested_workers,
            max_record_bytes,
            max_group_records,
            50,
        )
    }
    #[cfg(feature = "experimental-large")]
    pub fn select_experimental(
        budget: u64,
        requested_workers: u64,
        max_record_bytes: u64,
        max_group_records: u64,
    ) -> Result<(Self, Plan)> {
        if ![100, 250, 500].contains(&max_group_records) {
            return Err(invalid("admission: unsupported experimental group limit"));
        }
        Self::select_bounded(
            budget,
            requested_workers,
            max_record_bytes,
            max_group_records,
            500,
        )
    }
    fn select_bounded(
        budget: u64,
        requested_workers: u64,
        max_record_bytes: u64,
        max_group_records: u64,
        ceiling: u64,
    ) -> Result<(Self, Plan)> {
        if requested_workers == 0
            || requested_workers > 256
            || max_record_bytes == 0
            || max_record_bytes > PARSER_LINE
            || max_group_records == 0
            || max_group_records > ceiling
        {
            return Err(invalid("admission: invalid selection bounds"));
        }
        let min_chunk = add(mul(add(max_record_bytes, 16)?, 42)?, 4096)?.max(MIB);
        for workers in (1..=requested_workers).rev() {
            for (chunk, fan_in) in [(min_chunk.max(8 * MIB), 4), (min_chunk, 2)] {
                let request = Self {
                    budget,
                    workers,
                    window: workers,
                    merge_fan_in: fan_in,
                    max_record_bytes,
                    max_group_records,
                    sort_chunk_bytes: chunk,
                };
                if let Ok(plan) = request.plan_bounded(ceiling) {
                    return Ok((request, plan));
                }
            }
        }
        Err(invalid(
            "admission: budget below minimum selectable configuration",
        ))
    }
}
