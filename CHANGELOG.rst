
Changelog
=========

2.4.0
-----

- The ``--no-indels`` flags only discards indels, it cannot enable their
  analysis
- Changed the flag for signature correction
- Compress output with gzip. If output name is provided and does not exits,
  only the tsv file with the output is created.
- Added seed as option in the configuration file
- Updated dependencies and dropped Python 3.5 support

2.3.0
-----

- Added option to set the size limit for the indels

2.2.0
-----

- Simplified signature computation using bgsignature.
  Complex options should be computed externally.

- Adapted to bgparsers 0.9 that require a regions file with header.

- Simplified indels computation as assuming everything in the positive strand.

- Set by default the indels computation as substitutions.