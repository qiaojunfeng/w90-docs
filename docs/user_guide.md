---
author:
- Version 3.1
bibliography:
- ../wannier90.bib
title: "`wannier90`: User Guide"
---


[^1]: Technically, this is true for the case of an isolated group of $N$
    bands from which we obtain $N$ MLWF. When using the disentanglement
    procedure of Ref. [@souza-prb01], $\mathbf{A}^{(\mathbf{k})}$, for
    example, is a rectangular matrix. See
    Section [1.1](#sec:disentangle){reference-type="ref"
    reference="sec:disentangle"}.

[^2]: As ${\bf U}^{{\rm dis}({\bf k})}$ is a rectangular matrix this is
    a unitary operation in the sense that $({\bf U}^{{\rm
     dis}({\bf k})})^{\dagger}{\bf U}^{{\rm dis}({\bf k})}={\bf 1}_N$.

[^3]: It's worth noting that another visualisation program, VMD
    (<http://www.ks.uiuc.edu/Research/vmd>), is able to deal with
    certain special cases of non-orthogonal lattice vectors; see
    <http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html>
    for details.

[^4]: Scanning the Fermi level is currently supported only by the
    `postw90` module `berry`, for `berry_task=ahc,morb`. For all other
    functionalities that require a knowledge of $\varepsilon_F$, use
    `fermi_energy` instead.

[^5]: Note that there is a small bug with this feature in v3.2 (and
    subsequent patches) of ` quantum-espresso`. Please use a later
    version (if available) or the CVS version of `pw2wannier90.f90`,
    which has been fixed.

[^6]: Note that in `BoltzWann` the adaptive (energy) smearing scheme
    also implements a simple adaptive $k-$mesh scheme: if at any given
    $k$ point one of the band gradients is zero, then that $k$ point is
    replaced by 8 neighboring $k$ points. Thus, the final results for
    the DOS may be slightly different with respect to that given by the
    `dos` module.
