### [`create_coadd.py`](./create_coadd.py)

- `util.create_coadd1d_file`: use this in a reduction folder (like `all`) to create individual `*.coadd1d` file for each target.

```python
from highz_qso_arxiv.util import create_coadd1d_file
create_coadd1d_file("keck_lris_red_mark4.coadd1d", "GD153_lris_sens.fits")
```

- `util.create_coadd2d_file`:  use this in a reduction folder (like `all`) to create individual `*_coadd2d.cfg` file for each target in `coadd2d` folder.

```python
from highz_qso_arxiv.util import create_coadd2d_file
create_coadd2d_file("keck_mosfire_D.pypeit", spectrograph="keck_mosfire")
```

### [`util.py`](./util.py)