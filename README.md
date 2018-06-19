# BioSig
Genomic Sequence Indexing/Signature Generation and Searching.

A partial revive and simplification of Dr. Timothy Chappell's [TopSig](https://github.com/tachappell/topsig), focusing on genomic documents.

## Getting Started
Utilising the `Makefile`, running `make` will build an executable file `biosig` that can be used as follows.

## Indexing (`index`)
Generate a signature file from the given sequence file/s.

`./biosig index [OPTIONS] sequenceFile [sequenceFile2 [...]] -o outfile.bsig`

Available `OPTIONS`:

| Option      | Default    | Description                                                                                                      |
| ----------- | ---------- | ---------------------------------------------------------------------------------------------------------------- |
| -kmerlen    | 5          | Kmer length to hash.                                                                                             |
| -sigwidth   | 1024       | Signature size in bits.                                                                                          |
| -sigdensity | 19         | Signature density `1/x`.                                                                                         |
| -match      | (disabled) | Match and store sequence IDs with the given regular expression.<br/>Will prioritise first group `()` if present. |

Indexing will result in the output of both a signature file `outfile.bsig`, as well as a header file `outfile.bsig.head`, with the following formats:

**Signature File `.bsig`** (continuous binary)

This file contains the signatures that were generated from the sequences.
```
[signature A][signature B][...]
```

**Header File `.head`** (newline-delimited plaintext)

This file contains the IDs of the sequences that were indexed. The first line contains metadata regarding the indexing parameters.
```
[kmer length],[signature width],[signature density]
[sequence id A]
[sequence id B]
[...]
```

## Searching (`search`)
Search/Compare the signatures in a given signature file with one or more query signature files.

`./biosig search [OPTIONS] targetSig.bsig querySig.bsig [querySig2.bsig [...]] -o resultfile`

Available `OPTIONS`:

| Option      | Default      | Description                                                                                                      |
| ----------- | ------------ | ---------------------------------------------------------------------------------------------------------------- |
| -threshold  | 0.0          | Filter results lower than the given threshold `0-1`.                                                             |
| -top        | 0 (disabled) | Retain only the top k results.                                                                                   |
| -format     | tsv          | Result output format. `(tsv \| csv \| trec)`                                                                     |
| -unique     | (disabled)   | Do not compare signatures with identical IDs.                                                                    |
| -match      | (disabled)   | Match and store sequence IDs with the given regular expression.<br/>Will prioritise first group `()` if present. |

Results will be output as follows, depending on the format:

**tsv** (tab-separated values)
```
[query]    [target]    [hamming distance]    [normalised distance]
```

**csv** (comma-separated values)
```
[query],[target],[hamming distance],[normalised distance]
```

**trec** (Text REtrieval Conference)
```
[query] Q0 [target] [rank] [normalised distance] biosig
```
