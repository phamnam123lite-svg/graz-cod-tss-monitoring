# Repository layout

Use the root README for commands on main. All Python entry points are in scripts/. They resolve source and generated data relative to the repository root, not the current working directory.

Results in results/archives/ preserve their original internal paths and bytes. Extract them at the root to read analysis/ and supplemental_analysis/. Neither the data selection nor numerical methods changed in this restructuring.

verification/SHA256_manifest.json uses paths relative to the repository root. Older audit logs record the layout at the time they ran and should not be copied as current commands. Current commands are in README.md. Manuscript revision 5 cites release v1.3.0, which uses this layout and adds revision_analysis/ (see docs/REVISION_5.md). Tags v1.1.0 and v1.2.0 keep the earlier flat layout.
