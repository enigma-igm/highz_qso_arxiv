### MOSFIRE PROCEDURES

- Adjust `.pypeit` file, make use of `pypeit_view_fits` to check the dithering.
- Run `.pypeit` file.
- Generate sensitivity function with `pypeit_sensfunc`.
- Create coadd2d files with `highz_qso_arxiv.util.create_coadd2d_file`.
- Run all the coadd2d with shell command (for loop).
- Check the coadd2d files with `pypeit_show_2dspec`. If there is no 2d spectrum, there might be failed frames or no standard star. For the first case, find the failed frames and delete them. For the second case, provide offsets manually.
- If the target is not detected but we can see the trace on 2d, adjust the following arguments in the coadd2d configuration files (maybe we can change the `find_trim_edge` when generating the coadd2d files):

```
[reduce] 
    [[findobj]] 
        snr_thresh=1
        find_trim_edge = 600,600
        maxnumber_sci = 20
```

- Sometimes the names of the objects are not in Jxxxx+xxxx form (it would be better to change them into this form in `.pypeit`), and we cannot use `create_coadd1d_file`. In this case, we can use script to help us generate 1d files:

```python
def extract_from_coadd1d(target, seq):
  if type(seq) == int:
      seq = str(seq)
  outfile = open(f'{target}.coadd1d', 'w+')
  for idx, line in enumerate(lines):
      if idx < 14 or idx > 326:
          if idx == 5:
              i = line.find('_coadd')
              line = line[:i] + target + line[i:]
          print(line, end='', file=outfile)
      elif seq in line:
          print(line, end='', file=outfile)
  outfile.close()
```

- Check coadd2d configuration files in turn and make notes. Redo the coadd2d for those who need. After all the iterations finished, we do coadd1d.