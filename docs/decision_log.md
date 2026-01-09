# Decision Log (append-only)

## 2026-01-09 — Execution model: Container-as-a-Function
**Decision:** Tools in this repo must follow the Container-as-a-Function (CAF) model: non-interactive execution, explicit inputs via CLI/config, and outputs written only under a user-specified output directory (default `/output` in containers).
**Rationale:** Reproducible, HPC-friendly execution and consistent usage across tool repos.
**Consequences:** One canonical entrypoint; no interactive shells; no hidden state; stable output layout.

## 2026-01-09 — Artifact distribution policy
**Decision:** Do not commit container binaries to git. The canonical runtime artifact is a Docker/OCI image published to GHCR. Singularity/Apptainer is supported by pulling/converting from the Docker image (`singularity pull ... docker://...`). Optional prebuilt `.sif` files may be provided as GitHub Release assets with checksums.
**Rationale:** Avoid repo bloat; enable repeatable distribution; support HPC environments.
**Consequences:** Repos contain container specs + scripts + docs, not images. Release process may generate `.sif` artifacts.

## 2026-01-09 — Platform target for builds
**Decision:** Canonical build/target platform is `linux/amd64` to match Biowulf (x86_64).
**Rationale:** Reduce surprises with compiled dependencies and ensure compatibility with Biowulf.
**Consequences:** Prefer building containers on an x86_64 host (team VM). Apple Silicon builds must explicitly target amd64.

## 2026-01-09 — Documentation policy
**Decision:** `docs/session_notes.md` and `docs/decision_log.md` are committed and append-only. Session notes updated at least 1x/day during active work.
**Rationale:** Preserve context and prevent reintroducing rejected approaches.
**Consequences:** Major changes consult the decision log; sessions end with a notes update + commit/push prompt.
