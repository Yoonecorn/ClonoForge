#!/usr/bin/env python3
"""Independently verify a completed paired directory without original raw inputs."""
import argparse
import hashlib
import json
from pathlib import Path


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'duplicate JSON field: {key}')
        result[key] = value
    return result


def load(path):
    with path.open() as stream:
        return json.load(stream, object_pairs_hook=unique_object)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def fastq_count(path):
    count = 0
    with path.open('rb') as stream:
        while header := stream.readline():
            sequence, plus, quality = (stream.readline() for _ in range(3))
            if not all(x.endswith(b'\n') for x in (header, sequence, plus, quality)):
                raise ValueError(f'incomplete FASTQ: {path}')
            sequence, quality = sequence[:-1], quality[:-1]
            if (not header.startswith(b'@') or plus != b'+\n' or not sequence
                    or len(sequence) != len(quality)
                    or any(b not in b'ACGTN' for b in sequence)
                    or any(q < 33 or q > 126 for q in quality)):
                raise ValueError(f'invalid FASTQ: {path}')
            count += 1
    return count


def verify(root):
    receipt = load(root / 'completion.json')
    if receipt['schema'] not in ('clonoforge-m4-preliminary-v1', 'clonoforge-paired-v2') or receipt['complete'] is not True:
        raise ValueError('unsupported or incomplete receipt')
    expected_names = {'qualified.fastq', 'm3_completion.json', 'm3_consensus-pass.fastq'}
    if set(receipt['files']) != expected_names:
        raise ValueError('unexpected completion file set')
    for name, meta in receipt['files'].items():
        path = root / name
        if meta is None:
            if name != 'm3_consensus-pass.fastq' or path.exists() or path.is_symlink():
                raise ValueError(f'unexpected absent output: {name}')
            continue
        if (path.is_symlink() or not path.is_file() or path.stat().st_size != meta['bytes']
                or digest(path) != meta['sha256']):
            raise ValueError(f'file size/hash mismatch: {name}')
    counts = receipt['counts']
    if counts['qualified'] + counts['group_excluded'] != counts['candidates']:
        raise ValueError('candidate accounting mismatch')
    terminal = ['invalid_key', 'short_target', 'primer_miss', 'too_short', 'too_many_n', 'low_quality', 'candidates']
    if receipt['schema'] == 'clonoforge-paired-v2':
        terminal += ['assembly_no_overlap', 'assembly_ambiguous']
    if sum(counts[key] for key in terminal) != counts['pairs']:
        raise ValueError('pair accounting mismatch')
    if fastq_count(root / 'qualified.fastq') != counts['qualified']:
        raise ValueError('qualified FASTQ count mismatch')
    nested = load(root / 'm3_completion.json')
    if nested['success'] is not True:
        raise ValueError('nested M3 unsuccessful')
    total = 0
    for entry in nested['files']:
        if entry['name'] != 'm3_consensus-pass.fastq':
            raise ValueError('unexpected M3 category')
        path = root / entry['name']
        if entry['bytes'] != path.stat().st_size or entry['sha256'] != digest(path):
            raise ValueError('M3 metadata mismatch')
        actual = fastq_count(path)
        if actual != entry['records']:
            raise ValueError('consensus record count mismatch')
        total += actual
    if total != counts['qualified_groups']:
        raise ValueError('consensus group count mismatch')
    return {'status': 'PASS', 'completion_sha256': digest(root / 'completion.json'), 'counts': counts}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    print(json.dumps(verify(parser.parse_args().directory), sort_keys=True))
