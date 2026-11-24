# SDP Instance Space Analysis — Thesis Documents

LaTeX sources and notes for the thesis on **Instance Space Analysis (ISA) applied
to Semidefinite Programming (SDP)** problems.

The desktop application lives in its own repository:
[sdp-instance-space-thesis](https://github.com/Lucas-Zamora-F/sdp-instance-space-thesis).

## Structure

```
manual/      # User manual for the app
notes/       # Working notes (e.g. the CBLIB pipeline)
thesis 3/    # Thesis report 3
thesis 4/    # Thesis report 4
```

Build with LaTeX Workshop (recipe in `.vscode/settings.json`):
pdflatex → biber → makeglossaries → pdflatex × 2.

## License

See [LICENSE](LICENSE).
