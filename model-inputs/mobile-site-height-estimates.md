# Planning-derived antenna-height estimates

`mobile-site-height-estimates.csv` aligns one row to each of the 7,820 ComReg operator/site records, in the same order as `dist/data/mobile-sites.json`.

For 1,605 rows, an indexed planning description matched a ComReg record within 30 m and explicitly described a telecommunications support structure height. The model estimates the highest/top-mounted antenna height as that stated structure height minus 1 m. This is an estimate, not an antenna height measured in the field or confirmed as built. Rows retain planning status, authority, application number, match distance, source URL, quote and ComReg Site Identity for provenance.

The 1,605 estimates comprise 919 `existing_as_described`, 618 `proposed_or_replacement`, 44 `status_ambiguous` and 24 `planning_status_unspecified`. The proposed and ambiguous cases are included as a planning-informed scenario, not represented as current installed equipment. Most site/application identity links are automated spatial/text candidates; they have not each been manually reviewed.

The model retains the nominal 30 m transmitter-height default for the other 6,215 records. It excludes 110 alternate planning leads whose height estimate comes from a different point within 50 m, rather than the selected close match. No value is fabricated for a missing height.

The national planning index is published under CC BY 4.0: <https://data.gov.ie/dataset/national-planning-applications>. Local-authority descriptions are cited individually in the CSV when used. ComReg site coordinates and identities are from the Q1 2026 non-confidential mobile licence schedules: <https://www.comreg.ie/industry/radio-spectrum/licensing/search-licence-type/mobile-licences-2/>.

This inventory is not a complete or as-built mast-height census. The resulting signal map remains indicative and uncalibrated; actual antenna configuration, active EIRP, sector azimuth and downtilt remain unknown.
