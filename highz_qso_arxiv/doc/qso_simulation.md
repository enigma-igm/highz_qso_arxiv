Determine a signal-to-noise threshold for quasar confirmation (emperically?), then we use simulation tool to see how much exposure time we need to confirm a quasar with given magnitude and redshift. We can scale this to other magnitude with CCD equation.

Q: why can't we scale the $t_{\rm exp}$ with CCD equation?

Q: the sky noise is from a 300 sec real data. when we use larger exporsure time in simulation, we should also scale the real data?

<center>
<img src="./plot/sim.png" width="400"/>
</center>

CCD equation:

$$
\begin{aligned}
    \frac{S}{N}&=\frac{N_*}{\sqrt{N_*+n_{\rm pix}(N_S+N_D+N_R^2)}}\\
    &\approx \frac{N_*}{\sqrt{n_{\rm pix}N_S}}\\
    &\approx \frac{n_* t}{\sqrt{n_{\rm pix} n_S t}}\propto\sqrt{t}
\end{aligned}
$$

