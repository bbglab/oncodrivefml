
Changelog
=========

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