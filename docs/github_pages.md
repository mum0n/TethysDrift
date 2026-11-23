# Interactive Dashboard: Local Viewing and GitHub Pages Deployment


```

git switch gh-pages
Copy-Item ..\lakehouse\TethysDrift\snowcrab\interactive_larval_tracks.html .\index.html
Copy-Item ..\lakehouse\TethysDrift\snowcrab\interactive_larval_tracks.data.zz .\interactive_larval_tracks.data.zz
git add index.html interactive_larval_tracks.data.zz
git commit --amend --no-edit
git push --force-with-lease origin gh-pages
git switch main

```



`export_interactive_tracks_html` writes two paired files to the configured output directory (for example, `../lakehouse/TethysDrift/snowcrab/`):

- `interactive_larval_tracks.html`: Canvas-backed (`preferCanvas: true`) Leaflet.js dashboard shell with collapsible left controls, top-right palette/legend dock, single-click `n`-nearest release divergence inspection, and hover tooltips documenting controls and habitat assumptions.
- `interactive_larval_tracks.data.zz`: RFC 1950 zlib-compressed JSON sidecar containing daily-resolution particle trajectories (`max_points_per_track = 550` in `driver.jl`), developmental stages (including pre-hatch `:egg`), and settlement density.

Eulerian hydrodynamic grids are excluded from the dashboard payload (handled by `TethysFlow`). None of the `.zarr` archives or `.nc` files are needed by the browser.

## Opening Locally

- **Direct open (`file://`)**:
  - Open `interactive_larval_tracks.html` directly in a browser.
  - Because browsers block `fetch()` requests over `file://`, the loader overlay displays a **Select `interactive_larval_tracks.data.zz`** button.
  - Select `interactive_larval_tracks.data.zz` from the same directory; the browser decompresses the stream in memory via `DecompressionStream("deflate")` and renders the dashboard.
- **Local static HTTP server (automatic loading)**:
  - Start a static HTTP server in the directory containing both files:
    ```powershell
    cd ../lakehouse/TethysDrift/snowcrab
    python -m http.server 8000
    ```
  - Open `http://localhost:8000/interactive_larval_tracks.html`. The page fetches and decompresses `./interactive_larval_tracks.data.zz` automatically, matching GitHub Pages behavior.

## Publishing with an Orphan `gh-pages` Branch

An orphan `gh-pages` branch keeps binary `.data.zz` sidecars out of the `main` branch history.

- **Create the orphan branch and clear the tracked index** (from the `TethysDrift` repository root with a clean working tree on `main`):
  ```powershell
  git switch --orphan gh-pages
  git rm -rf .
  ```
  - `.git` and any files that were untracked or gitignored on `main` (such as `Manifest.toml` or local `test/` artifacts) remain on disk and should be left untouched so they persist when switching back to `main`.
- **Copy the dashboard files and add `.nojekyll`**:
  ```powershell
  Copy-Item ../lakehouse/TethysDrift/snowcrab/interactive_larval_tracks.html ./index.html
  Copy-Item ../lakehouse/TethysDrift/snowcrab/interactive_larval_tracks.data.zz ./interactive_larval_tracks.data.zz
  New-Item -ItemType File -Name .nojekyll
  ```
  - For multiple species or runs, place each pair in its own subdirectory (for example, `./snowcrab/interactive_larval_tracks.html` and `./snowcrab/interactive_larval_tracks.data.zz`).
- **Stage explicitly, commit, push, and return to `main`**:
  ```powershell
  git add index.html interactive_larval_tracks.data.zz .nojekyll
  git commit -m "Publish snow crab larval dispersal dashboard"
  git push -u origin gh-pages
  git switch main
  ```
  - Avoid `git add .` on the orphan branch, because `.gitignore` is absent there and untracked files such as `Manifest.toml` would be staged.
- **Enable GitHub Pages**:
  - In the repository settings on GitHub, open **Settings -> Pages -> Build and deployment**, set **Source** to **Deploy from a branch**, select **`gh-pages`** and **`/ (root)`**, and save.

## Updating Published Results

- Switch to `gh-pages`, copy the regenerated `.html` and `.data.zz` files, and stage them explicitly:
  ```powershell
  git switch gh-pages
  Copy-Item ../lakehouse/TethysDrift/snowcrab/interactive_larval_tracks.html ./index.html
  Copy-Item ../lakehouse/TethysDrift/snowcrab/interactive_larval_tracks.data.zz ./interactive_larval_tracks.data.zz
  git add index.html interactive_larval_tracks.data.zz
  ```
- Either append a normal commit or replace the previous snapshot in history to prevent repository growth:
  ```powershell
  git commit --amend --no-edit
  git push --force-with-lease origin gh-pages
  git switch main
  ```

## Alternative Layout (`docs/` on `main`)

- If maintaining a separate branch is unnecessary, place the paired `.html` and `.data.zz` files alongside an empty `docs/.nojekyll` under `docs/` (or a subdirectory such as `docs/snowcrab/`), commit on `main`, and configure GitHub Pages to deploy from **`main`** / **`/docs`**.

## Deployment Constraints

- Keep the `.html` and `.data.zz` files in the same directory and retain the original `.data.zz` filename (`interactive_larval_tracks.html` fetches `./interactive_larval_tracks.data.zz`; renaming the `.html` file to `index.html` is supported).
- Commit `.data.zz` with standard git rather than Git LFS, as GitHub Pages serves raw LFS pointer files rather than the binary payload.
- Individual files must remain below GitHub's 100 MB per-file limit.
