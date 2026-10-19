# Third-party software and notices

Acquisition, pins, configuration and evidence levels are in [EXTERNAL_DEPENDENCIES.md](EXTERNAL_DEPENDENCIES.md). [follow-up licensing decisions](docs/FOLLOW_UP_ITEMS.md) explains the distribution boundaries; no root license has been granted.

| Component | License evidence | Included in candidate? |
|---|---|---|
| MetaKSSD | Local upstream LICENSE.txt: Apache-2.0; source/hash recorded | No upstream source/binary; external query wrapper only |
| FastSpar v1.0.0 | Official versioned README: GPLv3 | No upstream source/binary; invocation and MACK-DB postprocessing only |
| fastCCLasso | Official GNU LGPL-3.txt: LGPLv3; exact core hash recorded | No official R core; wrapper loads explicitly configured user copy |
| Upstream DKI/cNODE ancestry | MIT, copyright 2024 Xu-Wen Wang | Adapted model code present; original notice retained at docs/third_party/DKI_UPSTREAM_MIT_LICENSE.txt |
| vegan / microeco / igraph | Current metadata: GPL-2 / GPL-3 / GPL>=2 | Installed packages called externally, no package source vendored |
| optparse / tidyverse / data.table | Current metadata: GPL>=2 / MIT+LICENSE / MPL-2.0 or LICENSE | External dependencies |
| PyTorch / torchdiffeq and other plotting/statistical packages | Existing environment records; obtain each package's official notices before redistribution | No package distribution bundled |

The installed metadata and archived versions have different evidence meanings; they are not interchangeable. Full local probe files are retained outside this Git checkout. Public upstream pages and local official checkout files were inspected read-only; no package install, upstream inference or AWS operation occurred. External-author research data also has separate rights and access conditions.
