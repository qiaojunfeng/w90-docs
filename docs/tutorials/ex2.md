# 35: Silicon — Projectability-disentangled Wannier functions with custom projectors

-   Outline: *Obtain MLWFs for silicon using projectability
    disentanglement, with additional $3d$ projectors to describe
    high-energy conduction bands. For more details on the methodology,
    see Ref. [@Qiao2023-pdwf]*

-   Directory: `examples/example35/`

-   Input Files

    -    `silicon.scf` *The  input file for ground state calculation*

    -    `silicon.bands` *The  input file for band structure
        calculation*

    -    `silicon.nscf` *The  input file to obtain Bloch states on a
        uniform grid*

    -    `silicon.pw2wan` *Input file for `pw2wannier90`*

    -    `silicon.win` *The `wannier90` input file*

1.  Run  to obtain the ground state of silicon\
    `pw.x < silicon.scf > scf.out`

2.  Run  to obtain the band structure of silicon\
    `pw.x < silicon.bands > bands.out`

3.  Run `bands.x` to obtain a `silicon.bands.dat` file containing the
    band structure of silicon\
    `bands.x < silicon.bandsx > bandsx.out`

4.  Run  to obtain the Bloch states on a uniform k-point grid\
    `pw.x < silicon.nscf > nscf.out`

5.  Run  to generate a list of the required overlaps (written into the
    `silicon.nnkp` file).\
    (note: see `win` input file, no need to specify initial projections,
    they are chosen from the pseudo-atomic orbitals inside the
    `ext_proj/Si.dat` file)\
    `wannier90.x -pp silicon`

6.  Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `silicon.mmn`
    and `silicon.amn` files).\
    `pw2wannier90.x < silicon.pw2wan > pw2wan.out`

7.  Run  to compute the MLWFs.\
    `wannier90.x silicon`

8.  Run `gnuplot` to compare DFT and Wannier-interpolated bands, this
    will generate a PDF file `silicon_bandsdiff.pdf`, see
    Fig. [1](#fig:silicon_bandsdiff){reference-type="ref"
    reference="fig:silicon_bandsdiff"}.\
    `./silicon_bandsdiff.gnu`

    <figure id="fig:silicon_bandsdiff">
    <div class="center">
    <embed src="silicon_bandsdiff.pdf" style="width:60.0%" />
    </div>
    <figcaption>Comparison of DFT and Wannier bands for
    silicon.</figcaption>
    </figure>

9.  (Optional) Clean up all output files\
    `make clean`

## Further ideas {#further-ideas .unnumbered}

1.  Try changing the `atom_proj_exclude` in `silicon.pw2wan` file, i.e.,
    these commented lines

    ` `

    > ! for excluding specific projectors\
    > ! this excludes 3d projectors, then the results are similar\
    > ! to that of using UPF file, i.e., project onto Si s+p orbitals\
    > ! for the indices of orbitals, see the pw2wan stdout\
    > ! atom_proj_exclude = 5 6 7 8 9 14 15 16 17 18\

2.  Now that $3d$ projectors provide us a larger space for optimization,
    you can try increasing the `dis_froz_max` to freeze higher energy
    bands, if you are targeting at reproducing those eigenvalues.

    Note that the `dis_proj_min/max` and `dis_froz_min/max` can be
    enabled simultaneously: the union of inner energy window and
    high-projectability states will be freezed, and the union of states
    outside outer energy window and having low projectability will be
    discared. Thus, you can still use energy window to make sure
    near-Fermi energy states are well reproduced, and use
    "projectability window" to selectively freeze atomic-like states in
    the conduction region.

3.  The default `dis_proj_max = 0.95` might not freeze all the states
    you want, try changing this value and see the band interpolation
    results. For other materials, it might worth trying decreasing this
    value to freeze more states.
