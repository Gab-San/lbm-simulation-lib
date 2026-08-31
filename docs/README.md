# Documentation

| Path | What it holds |
|------|---------------|
| `mainpage.md` | landing page of the Doxygen portal |
| `pages/` | the narrative pages: architecture, configuration, output formats, validation, performance, extending |
| `report/` | the project report: error results, profiling, flow animations, the UML sketch and the theory PDF |
| `report/assets/` | figures and error CSVs referenced by the report and by the pages |
| `imgResults/` | animations (GIF/MP4) of the 3D cases |
| `references/` | the Ghia et al. (1982) paper and the assignment proposal |
| `assets/` | the interactive D3Q19/D3Q27 direction viewers, copied verbatim into the HTML output |
| `CMakeLists.txt` | the `lbm-sim-docs` and `lbm-docs-full` Doxygen targets |

Build the portal with:

```bash
cmake --build build --target lbm-docs-full
```

Output lands in `docs/doxygen/`. The pages are written to be readable either
rendered by Doxygen or as plain Markdown in the repository.
