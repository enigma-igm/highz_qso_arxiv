### [`create_coadd1d_file.py`](./create_coadd1d_file.py)

- `util.create_coadd1d_file`: use this in a reduction folder (like `all`) to create individual `.coadd1d` file for each target.

```python
from highz_qso_arxiv.util import create_coadd1d_file
create_coadd1d_file("keck_lris_red_mark4.coadd1d", "GD153_lris_sens.fits")
```

### [`util.py`](./util.py)