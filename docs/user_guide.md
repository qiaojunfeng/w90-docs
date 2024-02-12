---
author:
- Version 3.1
bibliography:
- ../wannier90.bib
title: "`wannier90`: User Guide"
---

# Introduction

## Introduction {#introduction-1 .unnumbered}

### Getting Help {#getting-help .unnumbered}

The latest release of `wannier90` and documentation can always be found
at <http://www.wannier.org>.

The development version may be cloned/downloaded from the official
repository of the `wannier90` code on GitHub (see
<https://github.com/wannier-developers/wannier90>).

There is a `wannier90` mailing list for discussing issues in the
development, theory, coding and algorithms pertinent to MLWF. You can
register for this mailing list by following the links at
<http://www.wannier.org/forum.html>. Alternatively, for technical issues
about the `wannier90` code, check the official repository of
`wannier90` on GitHub where you may raise issues or ask questions about
about its functionalities.

Finally, many frequently asked questions are answered in
Appendix [2](#chap:faq){reference-type="ref" reference="chap:faq"}. An
expanded FAQ session may be found on the Wiki page of the GitHub
repository at
<https://github.com/wannier-developers/wannier90/wiki/FAQ>.

### Citation {#citation .unnumbered}

We ask that you acknowledge the use of `wannier90` in any publications
arising from the use of this code through the following reference

> \[ref\] G. Pizzi, V. Vitale, R. Arita, S. Blügel, F. Freimuth, G.
> Géranton, M. Gibertini, D. Gresch, C. Johnson, T. Koretsune, J.
> Ibañez-Azpiroz, H. Lee, J.M. Lihm, D. Marchand, A. Marrazzo, Y.
> Mokrousov, J.I. Mustafa, Y. Nohara, Y. Nomura, L. Paulatto, S. Poncé,
> T. Ponweiser, J. Qiao, F. Thöle, S.S. Tsirkin, M. Wierzbowska, N.
> Marzari, D. Vanderbilt, I. Souza, A.A. Mostofi, J.R. Yates,\
> Wannier90 as a community code: new features and applications, *J.
> Phys. Cond. Matt.* **32**, 165902 (2020)\
> <https://doi.org/10.1088/1361-648X/ab51ff>

If you are using versions 2.x of the code, cite instead:

> \[ref\] A. A. Mostofi, J. R. Yates, G. Pizzi, Y.-S. Lee, I. Souza,
> D. Vanderbilt and N. Marzari,\
> An updated version of `wannier90`: A Tool for Obtaining
> Maximally-Localised Wannier Functions, *Comput. Phys. Commun.*
> **185**, 2309 (2014)\
> <http://dx.doi.org/10.1016/j.cpc.2014.05.003>

It would also be appropriate to cite the original articles:\
\
Maximally localized generalized Wannier functions for composite energy
bands,\
N. Marzari and D. Vanderbilt, *Phys. Rev. B* **56**, 12847 (1997)\
\
Maximally localized Wannier functions for entangled energy bands,\
I. Souza, N. Marzari and D. Vanderbilt, *Phys. Rev. B* **65**, 035109
(2001)

### Credits {#credits .unnumbered}

The Wannier90 Developer Group includes Giovanni Pizzi (EPFL, CH),
Valerio Vitale (Cambridge, GB), David Vanderbilt (Rutgers University,
US), Nicola Marzari (EPFL, CH), Ivo Souza (Universidad del Pais Vasco,
ES), Arash A. Mostofi (Imperial College London, GB), and Jonathan R.
Yates (University of Oxford, GB).

The present release of `wannier90` was written by the Wannier90
Developer Group together with Ryotaro Arita (Riken and U. Tokyo, JP),
Stefan Blügel (FZ Jülich, DE), Frank Freimuth (FZ Jülich, DE), Guillame
Géranton (FZ Jülich, DE), Marco Gibertini (EPFL and University of
Geneva, CH), Dominik Gresch (ETHZ, CH), Charles Johnson (Imperial
College London, GB), Takashi Koretsune (Tohoku University and JST
PRESTO, JP), Julen Ibañez-Azpiroz (Universidad del Pais Vasco, ES),
Hyungjun Lee (EPFL, CH), Jae-Mo Lihm (Seoul National University, KR),
Daniel Marchand (EPFL, CH), Antimo Marrazzo (EPFL, CH), Yuriy Mokrousov
(FZ Jülich, DE), Jamal I. Mustafa (UC Berkeley, USA), Yoshiro Nohara
(Tokyo, JP), Yusuke Nomura (U. Tokyo, JP), Lorenzo Paulatto (Sorbonne
Paris, FR), Samuel Poncé (Oxford University, GB), Thomas Ponweiser (RISC
Software GmbH, AT), Florian Thöle (ETHZ, CH), Stepan Tsirkin
(Universidad del Pais Vasco, ES), Małgorzata Wierzbowska (Polish Academy
of Science, PL).

Contributors to the code include: Daniel Aberg (w90pov code), Lampros
Andrinopoulos (w90vdw code), Pablo Aguado Puente (gyrotropic routines),
Raffaello Bianco (k-slice plotting), Marco Buongiorno Nardelli (dosqc
v1.0 subroutines upon which transport.f90 is based), Stefano De
Gironcoli (pw2wannier90.x interface to Quantum ESPRESSO), Pablo Garcia
Fernandez (matrix elements of the position operator), Nicholas D. M.
Hine (w90vdw code), Young-Su Lee (specialised Gamma point routines and
transport), Antoine Levitt (preconditioning), Graham Lopez (extension of
pw2wannier90 to add terms needed for orbital magnetisation), Radu Miron
(constrained centres), Nicolas Poilvert (transport routines), Michel
Posternak (original plotting routines), Rei Sakuma (Symmetry-adapted
Wannier functions), Gabriele Sclauzero (disentanglement in spheres in
k-space), Matthew Shelley (transport routines), Christian Stieger
(routine to print the U matrices), David Strubbe (various
bugfixes/improvements), Timo Thonhauser (extension of pw2wannier90 to
add terms needed for orbital magnetisation), Junfeng Qiao (spin Hall
conductivity, projectability-disentangled Wannier functions).

We also acknowledge individuals not already mentioned above who
participated in the first Wannier90 community meeting (San Sebastian,
2016) for useful discussions: Daniel Fritsch, Victor Garcia Suarez,
Jan-Philipp Hanke, Ji Hoon Ryoo, Jürg Hutter, Javier Junquera, Liang
Liang, Michael Obermeyer, Gianluca Prandini, Paolo Umari.

`wannier90` Version 2.x was written by: Arash A. Mostofi, Giovanni
Pizzi, Ivo Souza, Jonathan R. Yates. `wannier90` Version 1.0 was written
by: Arash A. Mostofi, Jonathan R. Yates, Young-Su Lee. `wannier90` is
based on the Wannier Fortran 77 code written for isolated bands by
Nicola Marzari and David Vanderbilt and for entangled bands by Ivo
Souza, Nicola Marzari, and David Vanderbilt.

`wannier90` © 2007-2020 The Wannier Developer Group and individual
contributors

### Licence {#licence .unnumbered}

All the material in this distribution is free software; you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

# `wannier90.x`[]{#part:w90 label="part:w90"}

## Methodology {#sec:method}

`wannier90` computes maximally-localised Wannier functions (MLWF)
following the method of Marzari and Vanderbilt (MV) [@marzari-prb97].
For entangled energy bands, the method of Souza, Marzari and Vanderbilt
(SMV) [@souza-prb01] is used. We introduce briefly the methods and key
definitions here, but full details can be found in the original papers
and in Ref. [@mostofi-cpc08].

First-principles codes typically solve the electronic structure of
periodic materials in terms of Bloch states, $\psi_{n{\bf k}}$. These
extended states are characterised by a band index $n$ and crystal
momentum ${\bf k}$. An alternative representation can be given in terms
of spatially localised functions known as Wannier functions (WF). The WF
centred on a lattice site ${\bf R}$, $w_{n{\bf R}}({\bf r})$, is written
in terms of the set of Bloch states as

$$
w_{n{\bf R}}({\bf r})=\frac{V}{(2\pi)^3}\int_{\mathrm{BZ}}
\left[\sum_{m} U^{({\bf k})}_{mn} \psi_{m{\bf k}}({\bf
    r})\right]e^{-\mathrm{i}{\bf k}.{\bf R}} \:\mathrm{d}{\bf k} \ ,
$$

where $V$ is the unit cell volume, the integral is over the Brillouin
zone (BZ), and $\mathbf{U}^{(\mathbf{k})}$ is a unitary matrix that
mixes the Bloch states at each ${\bf k}$. $\mathbf{U}^{(\mathbf{k})}$ is
not uniquely defined and different choices will lead to WF with varying
spatial localisations. We define the spread $\Omega$ of the WF as
$$\Omega=\sum_n \left[\langle w_{n{\bf 0}}({\bf r})| r^2 | w_{n{\bf
      0}}({\bf r}) \rangle - | \langle w_{n{\bf 0}}({\bf r})| {\bf r}
      | w_{n{\bf 0}}({\bf r}) \rangle |^2 \right].$$ The total spread
can be decomposed into a gauge invariant term $\Omega_{\rm I}$ plus a
term ${\tilde \Omega}$ that is dependant on the gauge choice
$\mathbf{U}^{(\mathbf{k})}$. ${\tilde \Omega}$ can be further divided
into terms diagonal and off-diagonal in the WF basis, $\Omega_{\rm D}$
and $\Omega_{\rm OD}$,
$$\Omega=\Omega_{\rm I}+{\tilde \Omega}=\Omega_{\rm I}+\Omega_{\rm
  D}+\Omega_{\rm OD}$$ where
$$\Omega_{{\rm I}}=\sum_n \left[\langle w_{n{\bf 0}}({\bf r})| r^2 | w_{n{\bf
      0}}({\bf r}) \rangle - \sum_{{\bf R}m} \left| \langle w_{n{\bf
      R}}({\bf r})| {\bf r} | w_{n{\bf 0}}({\bf r}) \rangle \right| ^2
      \right]$$
$$\Omega_{\rm D}=\sum_n \sum_{{\bf R}\neq{\bf 0}} |\langle w_{n{\bf
    R}}({\bf r})| {\bf r} | w_{n{\bf 0}}({\bf r}) \rangle|^2$$
$$\Omega_{\rm OD}=\sum_{m\neq n} \sum_{{\bf R}} |\langle w_{m{\bf R}}({\bf
  r})| {\bf r} | w_{n{\bf 0}}({\bf r}) \rangle |^2$$ The MV method
minimises the gauge dependent spread $\tilde{\Omega}$ with respect the
set of $\mathbf{U}^{(\mathbf{k})}$ to obtain MLWF.

`wannier90` requires two ingredients from an initial electronic
structure calculation.

1.  The overlaps between the cell periodic part of the Bloch states
    $|u_{n{\bf k}}\rangle$
    $$M_{mn}^{(\bf{k,b})}=\langle u_{m{\bf k}}|u_{n{\bf k}+{\bf b}}\rangle,$$
    where the vectors ${\bf b}$, which connect a given k-point with its
    neighbours, are determined by `wannier90` according to the
    prescription outlined in Ref. [@marzari-prb97].

2.  As a starting guess the projection of the Bloch states
    $|\psi_{n\bf{k}}\rangle$ onto trial localised orbitals
    $|g_{n}\rangle$
    $$A_{mn}^{(\bf{k})}=\langle \psi_{m{\bf k}}|g_{n}\rangle,$$

Note that $\mathbf{M}^{(\mathbf{k},\mathbf{b})}$,
$\mathbf{A}^{(\mathbf{k})}$ and $\mathbf{U}^{(\mathbf{k})}$ are all
small, $N \times N$ matrices[^1] that are independent of the basis set
used to obtain the original Bloch states.

To date, `wannier90` has been used in combination with electronic codes
based on plane-waves and pseudopotentials (norm-conserving and
ultrasoft [@vanderbilt-prb90]) as well as mixed basis set techniques
such as FLAPW [@posternak-prb02].

### Entangled Energy Bands {#sec:disentangle}

The above description is sufficient to obtain MLWF for an isolated set
of bands, such as the valence states in an insulator. In order to obtain
MLWF for entangled energy bands we use the "disentanglement" procedure
introduced in Ref. [@souza-prb01].

We define an energy window (the "outer window"). At a given k-point
$\bf{k}$, $N^{({\bf k})}_{{\rm win}}$ states lie within this energy
window. We obtain a set of $N$ Bloch states by performing a unitary
transformation amongst the Bloch states which fall within the energy
window at each k-point:
$$| u_{n{\bf k}}^{{\rm opt}}\rangle = \sum_{m\in N^{({\bf k})}_{{\rm win}}}
U^{{\rm dis}({\bf k})}_{mn} | u_{m{\bf k}}\rangle$$ where
$\bf{U}^{{\rm dis}({\bf k})}$ is a rectangular
$N^{({\bf k})}_{{\rm win}} \times N$ matrix[^2]. The set of
$\bf{U}^{{\rm dis}({\bf k})}$ are obtained by minimising the gauge
invariant spread $\Omega_{{\rm I}}$ within the outer energy window. The
MV procedure can then be used to minimise $\tilde{\Omega}$ and hence
obtain MLWF for this optimal subspace.

It should be noted that the energy bands of this optimal subspace may
not correspond to any of the original energy bands (due to mixing
between states). In order to preserve exactly the properties of a system
in a given energy range (e.g., around the Fermi level) we introduce a
second energy window. States lying within this inner, or "frozen",
energy window are included unchanged in the optimal subspace.

## Parameters {#chap:parameters}

### Usage

`wannier90.x` can be run in parallel using MPI libraries to reduce the
computation time.

For serial execution use: `wannier90.x [-pp] [seedname]`

-    `seedname`: If a seedname string is given the code will read its
    input from a file `seedname.win`. The default value is `wannier`.
    One can also equivalently provide the string `seedname.win` instead
    of `seedname`.

-    `-pp`: This optional flag tells the code to generate a list of the
    required overlaps and then exit. This information is written to the
    file `seedname.nnkp`.

For parallel execution use:
`mpirun -np NUMPROCS wannier90.x [-pp] [seedname]`

-   `NUMPROCS`: substitute with the number of processors that you want
    to use.

Note that the `mpirun` command and command-line flags may be different
in your MPI implementation: read your MPI manual or ask your computer
administrator.

Note also that this requires that the `wannier90.x` executable has been
compiled in its parallel version (follow the instructions in the file
`README.install` in the main directory of the wannier90 distribution)
and that the MPI libraries and binaries are installed and correctly
configured on your machine.

### `seedname.win` File[]{#sec:seednamefile label="sec:seednamefile"}

The `wannier90` input file `seedname.win` has a flexible free-form
structure.

The ordering of the keywords is not significant. Case is ignored (so
`num_bands` is the same as `Num_Bands`). Characters after !, or \# are
treated as comments. Most keywords have a default value that is used
unless the keyword is given in `seedname.win`. Keywords can be set in
any of the following ways

` `

> num_wann 4
>
> num_wann = 4
>
> num_wann : 4

A logical keyword can be set to `true` using any of the following
strings: `T`, `true`, `.true.`.

For further examples see Section [10.1](#winfile){reference-type="ref"
reference="winfile"} and the the `wannier90` Tutorial.

### Keyword List {#parameter_data}

::: center
::: {#parameter_keywords1}
  ------------------- ------ ----------------------------------------------------------------------------------
        Keyword        Type  Description
                             
   System Parameters         
       num_wann         I    Number of WF
       num_bands        I    Number of bands passed to the code
    unit_cell_cart      P    Unit cell vectors in Cartesian coordinates
     atoms_cart \*      P    Positions of atoms in Cartesian coordinates
     atoms_frac \*      R    Positions of atoms in fractional coordinates with respect to the lattice vectors
        mp_grid         I    Dimensions of the Monkhorst-Pack grid of k-points
        kpoints         R    List of k-points in the Monkhorst-Pack grid
      gamma_only        L    Wavefunctions from underlying ab initio calculation are manifestly real
        spinors         L    WF are spinors
      shell_list        I    Which shells to use in finite difference formula
     search_shells      I    The number of shells to search when determining finite difference formula
     skip_B1_tests      L    Check the condition B1 of Ref. [@marzari-prb97]
        nnkpts          I    Explicit list of nearest-neighbour k-points.
       kmesh_tol        R    The tolerance to control if two kpoint belong to the same shell
  ------------------- ------ ----------------------------------------------------------------------------------

  : `seedname.win` file keywords defining the system. Argument types are
  represented by, I for a integer, R for a real number, P for a physical
  value, L for a logical value and S for a text string.\
  \* atoms_cart and atoms_frac may not both be defined in the same input
  file.
:::
:::

::: center
::: {#parameter_keywords2}
  ---------------------- ------ -------------------------------------------------------------------------------------------------
         Keyword          Type  Description
                                
       Job Control              
      postproc_setup       L    To output the `seedname.nnkp` file
      exclude_bands        I    List of bands to exclude from the calculation
    select_projections     I    List of projections to use in Wannierisation
     auto_projections      L    To automatically generate initial projections
         restart           S    Restart from checkpoint file
          iprint           I    Output verbosity level
       length_unit         S    System of units to output lengths
      wvfn_formatted       L    Read the wavefunctions from a (un)formatted file
           spin            S    Which spin channel to read
        devel_flag         S    Flag for development use
       timing_level        I    Determines amount of timing information written to output
       optimisation        I    Optimisation level
   translate_home_cell     L    To translate final Wannier centres to home unit cell when writing xyz file
        write_xyz          L    To write atomic positions and final centres in xyz file format
      write_vdw_data       L    To write data for futher processing by w90vdw utility
      write_hr_diag        L    To write the diagonal elements of the Hamiltonian in the Wannier basis to seedname.wout (in eV)
  ---------------------- ------ -------------------------------------------------------------------------------------------------

  : `seedname.win` file keywords defining job control. Argument types
  are represented by, I for a integer, R for a real number, P for a
  physical value, L for a logical value and S for a text string.
  translate_home_cell only relevant if write_xyz is `.true.`
:::
:::

::: center
::: {#parameter_keywords4}
  ---------------------------- ------ ---------------------------------------------------------------------------------------
            Keyword             Type  Description
                                      
   Disentanglement Parameters         
          dis_win_min            P    Bottom of the outer energy window
          dis_win_max            P    Top of the outer energy window
          dis_froz_min           P    Bottom of the inner (frozen) energy window
          dis_froz_max           P    Top of the inner (frozen) energy window
         dis_froz_proj           L    To activate projectability disentanglement
          dis_proj_min           P    Lower threshold for projectability disentanglement
          dis_proj_max           P    Upper threshold for projectability disentanglement
          dis_num_iter           I    Number of iterations for the minimisation of $\Omega_{\mathrm{I}}$
         dis_mix_ratio           R    Mixing ratio during the minimisation of $\Omega_{\mathrm{I}}$
          dis_conv_tol           R    The convergence tolerance for finding $\Omega_{\mathrm{I}}$
        dis_conv_window          I    The number of iterations over which convergence of $\Omega_{\mathrm{I}}$ is assessed.
        dis_spheres_num          I    Number of spheres in k-space where disentaglement is performed
     dis_spheres_first_wann      I    Index of the first band to be considered a Wannier function
          dis_spheres            R    List of centres and radii, for disentanglement only in spheres
  ---------------------------- ------ ---------------------------------------------------------------------------------------

  : `seedname.win` file keywords controlling the disentanglement.
  Argument types are represented by, I for a integer, R for a real
  number, P for a physical value, L for a logical value and S for a text
  string.
:::
:::

::: center
::: {#parameter_keywords5}
  ----------------------- ------ ------------------------------------------------------------------------------------------------------------------
          Keyword          Type  Description
                                 
   Wannierise Parameters         
         num_iter           I    Number of iterations for the minimisation of $\Omega$
       num_cg_steps         I    During the minimisation of $\Omega$ the number of Conjugate Gradient steps before resetting to Steepest Descents
        conv_window         I    The number of iterations over which convergence of $\Omega$ is assessed
         conv_tol           P    The convergence tolerance for finding $\Omega$
          precond           L    Use preconditioning
      conv_noise_amp        R    The amplitude of random noise applied towards end of minimisation procedure
      conv_noise_num        I    The number of times random noise is applied
      num_dump_cycles       I    Control frequency of check-pointing
     num_print_cycles       I    Control frequency of printing
        write_r2mn          L    Write matrix elements of $r^2$ between WF to file
      guiding_centres       L    Use guiding centres
     num_guide_cycles       I    Frequency of guiding centres
     num_no_guide_iter      I    The number of iterations after which guiding centres are used
       trial_step \*        R    The trial step length for the parabolic line search during the minimisation of $\Omega$
       fixed_step \*        R    The fixed step length to take during the minimisation of $\Omega$, instead of doing a parabolic line search
   use_bloch_phases \*\*    L    To use phases for initial projections
    site_symmetry\*\*\*     L    To construct symmetry-adapted Wannier functions
   symmetrize_eps\*\*\*     R    The convergence tolerance used in the symmetry-adapted mode
         slwf_num           I    The number of objective WFs for selective localization
      slwf_constrain        L    Whether to constrain the centres of the objective WFs
        slwf_lambda         R    Value of the Lagrange multiplier for constraining the objective WFs
       slwf_centres         P    The centres to which the objective WFs are to be constrained
  ----------------------- ------ ------------------------------------------------------------------------------------------------------------------

  : `seedname.win` file keywords controlling the wannierisation.
  Argument types are represented by, I for a integer, R for a real
  number, P for a physical value, L for a logical value and S for a text
  string. \* fixed_step and trial_step may not both be defined in the
  same input file. \*\*Cannot be used in conjunction with
  disentanglement. \*\*\*Cannot be used in conjunction with the inner
  (frozen) energy window.
:::
:::

::: {#parameter_keywords6}
  ------------------------------------------ ------ -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                   Keyword                    Type  Description
                                                    
               Plot Parameters                      
                 wannier_plot                  L    Plot the WF
              wannier_plot_list                I    List of WF to plot
            wannier_plot_supercell             I    Size of the supercell for plotting the WF
             wannier_plot_format               S    File format in which to plot the WF
              wannier_plot_mode                S    Mode in which to plot the WF, molecule or crystal
             wannier_plot_radius               R    Cut-off radius of WF\*
              wannier_plot_scale               R    Scaling parameter for cube files
           wannier_plot_spinor_mode            S    Quantity to plot for spinor WF
          wannier_plot_spinor_phase            L    Include the "phase" when plotting spinor WF
                  bands_plot                   L    Plot interpolated band structure
                 kpoint_path                   P    K-point path for the interpolated band structure
               bands_num_points                I    Number of points along the first section of the k-point path
              bands_plot_format                S    File format in which to plot the interpolated bands
              bands_plot_project               I    WF to project the band structure onto
               bands_plot_mode                 S    Slater-Koster type interpolation or Hamiltonian cut-off
                bands_plot_dim                 I    Dimension of the system
              fermi_surface_plot               L    Plot the Fermi surface
           fermi_surface_num_points            I    Number of points in the Fermi surface plot
                 fermi_energy                  P    The Fermi energy
               fermi_energy_min                P    Lower limit of the Fermi energy range
               fermi_energy_max                P    Upper limit of the Fermi energy range
              fermi_energy_step                R    Step for increasing the Fermi energy in the specified range
          fermi_surface_plot_format            S    File format for the Fermi surface plot
        [hr_plot]{style="color: red"}          L    [This parameter is not used anymore. Use write_hr instead.]{style="color: red"}
       [write_hr]{style="color: blue"}         L    [Write the Hamiltonian in the WF basis]{style="color: blue"}
      [write_rmn ]{style="color: blue"}        L    [Write the position operator in the WF basis]{style="color: blue"}
      [write_bvec ]{style="color: blue"}       L    [Write to file the matrix elements of the bvectors and their weights]{style="color: blue"}
       [write_tb ]{style="color: blue"}        L    [Write lattice vectors, Hamiltonian, and position operator in WF basis]{style="color: blue"}
                  hr_cutoff                    P    Cut-off for the absolute value of the Hamiltonian
                 dist_cutoff                   P    Cut-off for the distance between WF
               dist_cutoff_mode                S    Dimension in which the distance between WF is calculated
           translation_centre_frac             R    Centre of the unit cell to which final WF are translated
   [use_ws_distance ]{style="color: blue"}     L    [Improve interpolation using minimum distance between WFs, see Chap. [\[chap:interpolation\]](#chap:interpolation){reference-type="ref" reference="chap:interpolation"}]{style="color: blue"}
   [ws_distance_tol ]{style="color: blue"}     R    [Absolute tolerance for the distance to equivalent positions.]{style="color: blue"}
    [ws_search_size ]{style="color: blue"}     I    [Maximum extension in each direction of the super-cell of the Born-von Karmann cell to search for points inside the Wigner-Seitz cell]{style="color: blue"}
   [write_u_matrices ]{style="color: blue"}    L    [Write $\mathbf{U}^{(\mathbf{k})}$ and $\mathbf{U}^{\mathrm{dis}(\mathbf{k})}$ matrices to files]{style="color: blue"}
                                                    
  ------------------------------------------ ------ -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  : `seedname.win` file keywords controlling the plotting. Argument
  types are represented by, I for a integer, R for a real number, P for
  a physical value, L for a logical value and S for a text string. \*
  Only applies when wannier_plot_format is `cube`.
:::

::: center
::: {#parameter_keywords7}
  -------------------------- ------ -----------------------------------------------------------
           Keyword            Type  Description
                                    
     Transport Parameters           
          transport            L    Calculate quantum conductance and density of states
        transport_mode         S    Bulk or left-lead_conductor_right-lead calculation
         tran_win_min          P    Bottom of the energy window for transport calculation
         tran_win_max          P    Top of the energy window for transport calculation
       tran_energy_step        R    Sampling interval of the energy values
         fermi_energy          R    The Fermi energy
         tran_num_bb           I    Size of a bulk Hamiltonian
         tran_num_ll           I    Size of a left-lead Hamiltonian
         tran_num_rr           I    Size of a right-lead Hamiltonian
         tran_num_cc           I    Size of a conductor Hamiltonian
         tran_num_lc           I    Number of columns in a left-lead_conductor Hamiltonian
         tran_num_cr           I    Number of rows in a conductor_right-lead Hamiltonian
       tran_num_cell_ll        I    Number of unit cells in PL of left lead
       tran_num_cell_rr        I    Number of unit cells in PL of right lead
        tran_num_bandc         I    Half-bandwidth+1 of a band-diagonal conductor Hamiltonian
        tran_write_ht          L    Write the Hamiltonian for transport calculation
         tran_read_ht          L    Read the Hamiltonian for transport calculation
      tran_use_same_lead       L    Left and right leads are the same
     tran_group_threshold      R    Distance that determines the grouping of WFs
          hr_cutoff            P    Cut-off for the absolute value of the Hamiltonian
         dist_cutoff           P    Cut-off for the distance between WF
       dist_cutoff_mode        S    Dimension in which the distance between WF is calculated
         one_dim_axis          S    Extended direction for a one-dimensional system
   translation_centre_frac     R    Centre of the unit cell to which final WF are translated
  -------------------------- ------ -----------------------------------------------------------

  : `seedname.win` file keywords controlling transport. Argument types
  are represented by, I for a integer, R for a real number, P for a
  physical value, L for a logical value and S for a text string.
:::
:::

### System

#### `integer :: num_wann`

Number of WF to be found.

No default.

#### `integer :: num_bands`

Total number of bands passed to the code in the `seedname.mmn` file.

Default `num_bands`=`num_wann`

#### Cell Lattice Vectors

The cell lattice vectors should be specified in Cartesian coordinates.

`begin unit_cell_cart`\
`[units]` $$\begin{array}{ccc}
A_{1x} & A_{1y} & A_{1z} \\
A_{2x} & A_{2y} & A_{2z} \\
A_{3x} & A_{3y} & A_{3z}
\end{array}$$ `end unit_cell_cart`

Here $A_{1x}$ is the $x$-component of the first lattice vector
$\mathbf{A}_1$, $A_{2y}$ is the $y$-component of the second lattice
vector $\mathbf{A}_2$, etc.

`[units]` specifies the units in which the lattice vectors are defined:
either `Bohr` or `Ang`.

The default value is `Ang`.

#### Ionic Positions

The ionic positions may be specified in fractional coordinates relative
to the lattice vectors of the unit cell, or in absolute Cartesian
coordinates. Only one of `atoms_cart` and `atoms_frac` may be given in
the input file.

##### Cartesian coordinates

`begin atoms_cart`\
`[units]` $$\begin{array}{cccc}
P  & R^{P}_{x} & R^{P}_{y} & R^{P}_{z} \\
Q  & R^{Q}_{x} & R^{Q}_{y} & R^{Q}_{z} \\
\vdots
\end{array}$$ `end atoms_cart`

The first entry on a line is the atomic symbol. The next three entries
are the atom's position $\mathbf{R}=(R_x , R_y, R_z)$ in Cartesian
coordinates. The first line of the block, `[units]`, specifies the units
in which the coordinates are given and can be either `bohr` or `ang`. If
not present, the default is `ang`.

##### Fractional coordinates

`begin atoms_frac` $$\begin{array}{cccc}
P  & F^{P}_{1} & F^{P}_{2} & F^{P}_{3} \\
Q  & F^{Q}_{1} & F^{Q}_{2} & F^{Q}_{3} \\
\vdots
\end{array}$$ `end atoms_frac`

The first entry on a line is the atomic symbol. The next three entries
are the atom's position in fractional coordinates $\mathbf{F} = F_1
\mathbf{A}_{1} + F_2 \mathbf{A}_{2} + F_3 \mathbf{A}_{3}$ relative to
the cell lattice vectors $\mathbf{A}_i$, $i\in [1,3]$.

#### `integer, dimension :: mp_grid(3)`

Dimensions of the regular (Monkhorst-Pack) k-point mesh. For example,
for a $2\times2\times2$ grid:

`mp_grid : 2  2  2`

No default.

#### K-points

Each line gives the coordinate $\mathbf{K}=K_1 \mathbf{B}_{1} + K_2
\mathbf{B}_{2} + K_3 \mathbf{B}_3$ of a k-point in relative
(crystallographic) units, i.e., in fractional units with respect to the
primitive reciprocal lattice vectors $\mathbf{B}_{i}$, $i \in [1,3]$.
The position of each k-point in this list assigns its numbering; the
first k-point is k-point 1, the second is k-point 2, and so on.

`begin kpoints`\
$$\begin{array}{ccc}
 K^{1}_{1} & K^{1}_{2} & K^{1}_{3} \\
 K^{2}_{1} & K^{2}_{2} & K^{2}_{3} \\
\vdots
\end{array}$$ `end kpoints`

There is no default.

**Note**: There is an utility provided with `wannier90`, called
`kmesh.pl`, which helps to generate the explicit list of $k$ points
required by `wannier90`. See Sec. [1.1](#sec:kmesh){reference-type="ref"
reference="sec:kmesh"}.

#### `logical :: gamma_only`

If `gamma_only=true`, then `wannier90` uses a branch of algorithms for
disentanglement and localisation that exploit the fact that the Bloch
eigenstates obtained from the underlying ab initio calculation are
manifestly real. This can be the case when only the $\Gamma$-point is
used to sample the Brillouin zone. The localisation procedure that is
used in the $\Gamma$-only branch is based on the method of
Ref. [@gygi-cpc03].

The default value is `false`.

#### `logical :: spinors`

If `spinors=true`, then `wannier90` assumes that the WF correspond to
singularly occupied spinor states and `num_elec_per_state=1`.

The default value is `false`.

#### Shells

The MV scheme requires a finite difference expression for
$\nabla_{\bf k}$ defined on a uniform Monkhorst-Pack mesh of k-points.
The vectors $\{{\bf b}\}$ connect each mesh-point ${\bf k}$ to its
nearest neighbours. $N_{\mathrm{sh}}$ shells of neighbours are included
in the finite-difference formula, with $M_s$ vectors in the
$s^{\mathrm{th}}$ shell. For $\nabla_{{\bf k}}$ to be correct to linear
order, we require that the following equation is satisfied (Eq. B1 of
Ref. [@marzari-prb97]): $$\label{eq:B1}
\sum_{s}^{N_{\mathrm{sh}}} w_s \sum_i^{M_{\mathrm{s}}}
b_{\alpha}^{i,s} b_{\beta}^{i,s} = \delta_{\alpha\beta}\:,$$ where
${\bf b}^{i,s}$, $i\in[1,M_s]$, is the $i^{\mathrm{th}}$ vector
belonging to the $s^{\mathrm{th}}$ shell with associated weight $w_s$,
and $\alpha$ and $\beta$ run over the three Cartesian indices.

#### `integer :: shell_list(:)`

`shell_list` is vector listing the shells to include in the finite
difference expression. If this keyword is absent, the shells are chosen
automatically.

#### `integer :: search_shells`

Specifies the number of shells of neighbours over which to search in
attempting to determine an automatic solution to the B1 condition
Eq. [\[eq:B1\]](#eq:B1){reference-type="ref" reference="eq:B1"}. Larger
values than the default may be required in special cases e.g. for very
long thin unit cells.

The default value is 36.

#### `logical :: skip_B1_tests`

If set to `.true.`, does not check the B1 condition
Eq. [\[eq:B1\]](#eq:B1){reference-type="ref" reference="eq:B1"}. This
should *only* be used if one knows why the B1 condition should not be
verified. A typical use of this flag is in conjunction with the Z2PACK
code: <http://www.physics.rutgers.edu/z2pack/>.

The default value is `.false.`.

#### `integer, dimension(:, 5) :: nnkpts`

Specifies the nearest-neighbour k-points which are written to the
`.nnkp` file. This can be used to explicitly specify which overlap
matrices should be calculated.

    begin nnkpts
    1   2   0  0  0
    .
    .
    end nnkpts

Each nearest neighbour $\mathbf{k + b}$ is given by a line of 5
integers. The first specifies the k-point number `nkp` of $\mathbf{k}$.
The second is the k-point number of the neighbour. The final three
integers specify the reciprocal lattice vector which brings the k-point
specified by the second integer to $\mathbf{k + b}$.

This format is the same as in the `.nnkp` file, except that the number
of neighbours per k-point is not specified. However, the number of
neighbours still needs to be a multiple of the number of k-points.

This input parameter can be used only if `postproc_setup = .true.`, and
is not intended to be used with a full Wannier90 run. It can be used
also if the k-points do not describe a regular mesh.

#### `real(kind=dp) :: kmesh_tol`

Two kpoints belong to the same shell if the distance between them is
less than `kmesh_tol`. Units are Ang.

The default value is 0.000001 Ang.

### Projection

The projections block defines a set of localised functions used to
generate an initial guess for the unitary transformations. This data
will be written in the `seedname.nnkp` file to be used by a
first-principles code.

`begin projections`\
.\
.\
`end projections`

If `guiding_centres`=`true`, then the projection centres are used as the
guiding centres in the Wannierisation routine.

For details see Section [3.1](#sec:proj){reference-type="ref"
reference="sec:proj"}.

### Job Control

#### `logical :: postproc_setup`

If `postproc_setup`=`true`, then the wannier code will write
`seedname.nnkp` file and exit. If `wannier90` is called with the option
`-pp`, then `postproc_setup` is set to `true`, over-riding its value in
the `seedname.win` file.

The default value is `false`.

#### `integer :: iprint`

This indicates the level of verbosity of the output from 0 ("low"), the
bare minimum, to 3 ("high"), which corresponds to full debugging output.

The default value is 1.

#### `integer :: optimisation`

This indicates the level of optimisation used in the code. This is a
trade between speed and memory. A positive number indicates fastest
execution time at the cost of more memory. Zero or negative numbers
indicates a smaller memory footprint - at increased execution time.

At the moment the only values that have an effect are `optimisation<=0`
(low memory) and `optimisation>0` (fast)

The default value is 3.

#### `character(len=20) :: length_unit`

The length unit to be used for writing quantities in the output file
`seedname.wout`.

The valid options for this parameter are:

-   `Ang` (default)

-   `Bohr`

#### `character(len=50) :: devel_flag`

Not a regular keyword. Its purpose is to allow a developer to pass a
string into the code to be used inside a new routine as it is developed.

No default.

#### `integer :: exclude_bands(:)`

A k-point independent list of states to excluded from the calculation of
the overlap matrices; for example to select only valence states, or
ignore semi-core states. This keyword is passed to the first-principles
code via the `seedname.nnkp` file. For example, to exclude bands 2, 6,
7, 8 and 12:

`exclude_bands : 2, 6-8, 12`

#### `integer :: select_projections(:)`

A list of projections to be included in the wannierisation procedure. In
the case that `num_proj` is greater than `num_wann`, this keyword allows
a subset of the projections in the projection matrices to be used. For
example, to select the projections given by the indices 2, 6, 7, 8 and
12:

`select_projections : 2, 6-8, 12`

#### `logical :: auto_projections`

If `.true.` and no projections block is defined, then `wannier90` writes
an additional block in the `.nnkp` file during the pre-processing step,
to instruct the interface code to automatically generate the
$A_{mn}^{(\mathbf{k})}$.

For additional information on the behavior and on the added block, see
Sec. [\[sec:auto-projections-block\]](#sec:auto-projections-block){reference-type="ref"
reference="sec:auto-projections-block"}.

**Note:** the interface code (e.g. `pw2wannier90.x`) must have at least
one implementation of a method to automatically generate initial
projections in order for this option to be usable.

The default value of this parameter is $\verb#false#$.

#### `character(len=20) :: restart`

If `restart` is present the code will attempt to restart the calculation
from the `seedname.chk ` file. The value of the parameter determines the
position of the restart

The valid options for this parameter are:

-   `default`. Restart from the point at which the check file
    `seedname.chk` was written

-   `wannierise`. Restart from the beginning of the wannierise routine

-   `plot`. Go directly to the plotting phase

-   `transport`. Go directly to the transport routines

#### `character(len=20) :: wvfn_formatted`

If `wvfn_formatted`=`true`, then the wavefunctions will be read from
disk as formatted (ie ASCII) files; otherwise they will be read as
unformatted files. Unformatted is generally preferable as the files will
take less disk space and I/O is significantly faster. However such files
will not be transferable between all machine architectures and formatted
files should be used if transferability is required (i.e., for test
cases).

The default value of this parameter is $\verb#false#$.

#### `character(len=20) :: spin`

For bands from a spin polarised calculation `spin` determines which set
of bands to read in, either `up` or `down`.

The default value of this parameter is `up`.

#### `integer :: timing_level`

Determines the amount of timing information regarding the calculation
that will be written to the output file. A value of 1 produces the least
information.

The default value is 1.

#### `logical :: translate_home_cell`

Determines whether to translate the final Wannier centres to the home
unit cell at the end of the calculation. Mainly useful for molecular
systems in which the molecule resides entirely within the home unit cell
and user wants to write an xyz file (`write_xyz=.true.`) for the WF
centres to compare with the structure.

The default value is `false`.

#### `logical :: write_xyz`

Determines whether to write the atomic positions and final Wannier
centres to an xyz file, `seedname_centres.xyz`, for subsequent
visualisation.

The default value is `false`.

#### `logical :: write_vdw_data`

Determines whether to write `seedname.vdw` for subsequent
post-processing by the `w90vdw` utility (in the `utility/w90vdw/`
directory of the distribution) for calculating van der Waals energies.
Brillouin zone sampling must be at the Gamma-point only.

The default value is `false`.

### Disentanglement

These keywords control the disentanglement routine of
Ref. [@souza-prb01], i.e., the iterative minimisation of
$\Omega_{\mathrm{I}}$. This routine will be activated if
`num_wann`$\:<\:$`num_bands`.

#### `real(kind=dp) :: dis_win_min`

The lower bound of the outer energy window for the disentanglement
procedure. Units are eV.

The default is the lowest eigenvalue in the system.

#### `real(kind=dp) :: dis_win_max`

The upper bound of the outer energy window for the disentanglement
procedure. Units are eV.

The default is the highest eigenvalue in the given states (i.e., all
states are included in the disentanglement procedure).

#### `real(kind=dp) :: dis_froz_min`

The lower bound of the inner energy window for the disentanglement
procedure. Units are eV.

If `dis_froz_max` is given, then the default for `dis_froz_min` is
`dis_win_min`.

#### `real(kind=dp) :: dis_froz_max`

The upper bound of the inner (frozen) energy window for the
disentanglement procedure. If `dis_froz_max` is not specified, then
there are no frozen states. Units are eV.

No default.

#### `logical :: dis_froz_proj`

To activate projectability disentanglement procedure, which selectively
discard/disentangle/freeze state $\vert n \mathbf{k}\rangle$ based on
its projectability onto some localized atomic orbitals.

Note: this requires the `amn` file is properly normalized, i.e.,
projectability computed from $A A^\dagger$ must be smaller than or equal
to 1. The pseudo-atomic projection satisfies such requirement, see
[3.6](#sec:proj_pdwf){reference-type="ref" reference="sec:proj_pdwf"}.

Additionally, one can combine projectability disentanglement with energy
disentanglement, i.e., enable both `dis_proj_min/max` and
`dis_froz_min/max` simultaneously in the `win` file. These settings will
freeze the union of inner energy window and high-projectability states,
and exclude the union of states outside outer energy window and having
low projectability.

#### `real(kind=dp) :: dis_proj_min`

The lower bound for the projectability disentanglement procedure.

For states with projectabilities smaller than `dis_proj_min`, they will
be discarded in the disentanglement procedure, i.e., similar to the case
of outside of the outer energy window.

For states with projectabilities larger than or equal to `dis_proj_min`,
they will be included in the disentanglement procedure, i.e., similar to
the case of inside the outer energy window.

No unit.

The default value is 0.95.

#### `real(kind=dp) :: dis_proj_max`

The upper bound for the projectability disentanglement procedure. For
states with projectability larger than or equal to `dis_proj_max`, they
will be freezed in the disentanglement procedure, i.e., similar to the
case of inside the inner energy window.

No unit.

The default value is 0.01.

#### `integer :: dis_num_iter`

In the disentanglement procedure, the number of iterations used to
extract the most connected subspace.

The default value is 200.

#### `real(kind=dp) :: dis_mix_ratio`

In the disentanglement procedure, the mixing parameter to use for
convergence (see pages 4-5 of Ref. [@souza-prb01]). A value of 0.5 is a
'safe' choice. Using 1.0 (i.e., no mixing) often gives faster
convergence, but may cause the minimisation of $\Omega_{\mathrm{I}}$ to
be unstable in some cases.

Restriction: $0.0<\:$`dis_mix_ratio`$\:\leq 1.0$

The default value is 0.5

#### `real(kind=dp) :: dis_conv_tol`

In the disentanglement procedure, the minimisation of
$\Omega_{\mathrm{I}}$ is said to be converged if the fractional change
in the gauge-invariant spread between successive iterations is less than
`dis_conv_tol` for `dis_conv_window` iterations. Units are Å$^2$.

The default value is 1.0E-10

#### `integer :: dis_conv_window`

In the disentanglement procedure, the minimisation is said to be
converged if the fractional change in the spread between successive
iterations is less than `dis_conv_tol` for `dis_conv_window` iterations.

The default value of this parameter is 3.

#### `integer :: dis_spheres_num`

Number of spheres in reciprocal space where the k-dependent
disentanglement is performed. No disentanglement is performed for those
k-points that are not included in any of the spheres.

The default is 0, which means disentangle at every k-point in the full
BZ (the standard mode in Wannier90).

#### `integer :: dis_spheres_first_wann`

Index of the first band that has to be considered as a Wannier function.
Used only if `dis_spheres_num` is greater than zero. At k-points where
disentanglement is not performed the bands from `dis_spheres_first_wann`
to `dis_spheres_first_wann+num_wann` are used to wannierise. The bands
excluded using `exclude_bands` should not be counted.

The default is 1, the band at the lowest energy.

#### dis_spheres

Each line gives the coordinate $\mathbf{K}=K_1 \mathbf{B}_{1} + K_2
\mathbf{B}_{2} + K_3 \mathbf{B}_3$ of a k-point representing the center
of one of the spheres used for k-dependent disentanglement. The same
crystallographic units as for `kpoints` are used here. Each k-point
coordinate $\mathbf{K}^i$ must the followed by the respectice sphere
radius $r_{i}$ in inverse angstrom (on the same line).

The number of lines must be equal to `dis_spheres_num`.

`begin dis_spheres` $$\begin{array}{cccc}
 K^{1}_{1} & K^{1}_{2} & K^{1}_{3} & r_{1} \\
 K^{2}_{1} & K^{2}_{2} & K^{2}_{3} & r_{2} \\
\vdots
\end{array}$$ `end dis_spheres`

There is no default.

### Wannierise {#sec:wann_params}

Iterative minimisation of $\widetilde{\Omega}$, the non-gauge-invariant
part of the spread functional.

#### `integer :: num_iter`

Total number of iterations in the minimisation procedure. Set
`num_iter=0` if you wish to generate projected WFs rather than
maximally-localized WFs (see Example 8 in the Tutorial).

The default value is 100

#### `integer :: num_cg_steps`

Number of conjugate gradient steps to take before resetting to steepest
descents.

The default value is 5

#### `integer :: conv_window`

If `conv_window`$\:>1$, then the minimisation is said to be converged if
the change in $\Omega$ over ` conv_window` successive iterations is less
than ` conv_tol`. Otherwise, the minimisation proceeds for num_iter
iterations (default).

The default value is -1

#### `real(kind=dp) :: conv_tol`

If `conv_window`$\:>1$, then this is the convergence tolerance on
$\Omega$, otherwise not used. Units are Å$^2$.

The default value is 1.0E-10

#### `logical :: precond`

Whether or not to use preconditioning to speed up the minimization of
the spreads. This is based on the same idea as the classical
Tetter-Payne-Allan preconditionning for DFT and dampens the
high-frequency oscillations of the gradient due to contributions from
large real lattice vectors. It is useful when the optimization is slow,
especially on fine grids. When `optimisation<3`, this uses a slower
algorithm to save memory.

The default value is `false`.

#### `real(kind=dp) :: conv_noise_amp`

If `conv_noise_amp`$\:>0$, once convergence (as defined above) is
achieved, some random noise $f$ is added to the search direction, and
the minimisation is continued until convergence is achieved once more.
If the same value of $\Omega$ as before is arrived at, then the
calculation is considered to be converged. If not, then random noise is
added again and the procedure repeated up to a maximum of
` conv_noise_num` times. `conv_noise_amp` is the amplitude of the random
noise $f$ that is added to the search direction:
$0 < |f| <\:$`conv_noise_amp`. This functionality requires
` conv_window`$\:>1$. If `conv_window` is not specified, it is set to
the value 5 by default.

If `conv_noise_amp`$\:\leq 0$, then no noise is added (default).

The default value is -1.0

#### `integer :: conv_noise_num`

If `conv_noise_amp`$\:>0$, then this is the number of times in the
minimisation that random noise is added.

The default value is 3

#### `integer :: num_dump_cycles`

Write sufficient information to do a restart every `num_dump_cycles`
iterations.

The default is 100

#### `integer :: num_print_cycles`

Write data to the master output file `seedname.wout` every
`num_print_cycles` iterations.

The default is 1

#### `logical :: write_r2mn`

If $\verb#write_r2mn#=\verb#true#$, then the matrix elements
$\langle m|r^2|n\rangle$ (where $m$ and $n$ refer to WF) are written to
file `seedname.r2mn` at the end of the Wannierisation procedure.

The default value of this parameter is `false`.

#### `logical :: guiding_centres`

Use guiding centres during the minimisation, in order to avoid local
minima.

`wannier90` uses a logarithm definition of the spread functional. As we
are taking the log of a complex argument there is a possibility that the
algorithm might make inconsistent choices for the branch cut. This
manifests itself as complex WF with a large spread. By using guiding
centres the code will attempt to make a consistent choice of branch cut.
Experience shows that with `guiding_centres` set to true this problem is
avoided and doing so does not cause any problems. For this reason we
recommend setting `guiding_centres` to true where possible (it is only
not possible if an explicit projection block is not defined).

The default value is `false`.

#### `integer :: num_guide_cycles`

If `guiding_centres` is set to true, then the guiding centres are used
only every `num_guide_cycles`.

The default value is 1.

#### `integer :: num_no_guide_iter`

If `guiding_centres` is set to true, then the guiding centres are used
only after `num_no_guide_iter` minimisation iterations have been
completed.

The default value is 0.

#### `real(kind=dp) :: trial_step`

The value of the trial step for the parabolic fit in the line search
minimisation used in the minimisation of the spread function. Cannot be
used in conjunction with `fixed_step` (see below). If the minimisation
procedure doesn't converge, try decreasing the value of `trial_step` to
give a more accurate line search.

The default value is 2.0

#### `real(kind=dp) :: fixed_step`

If this is given a value in the input file, then a fixed step of length
`fixed_step` (instead of a parabolic line search) is used at each
iteration of the spread function minimisation. Cannot be used in
conjunction with `trial_step`. This can be useful in cases in which
minimisation with a line search fails to converge.

There is no default value.

#### `logical :: use_bloch_phases`

Determines whether to use the Bloch functions as the initial guess for
the projections. Can only be used if `disentanglement = false`.

The default value is `false`.

#### `logical :: site_symmetry`

Construct symmetry-adapted Wannier functions. For the detail of the
theoretical background, see Ref. [@sakuma-prb13]. Cannot be used in
conjunction with the inner (frozen) energy window.

The default value is `false`.

#### `real(kind=dp) :: symmetrize_eps`

Convergence threshold to check whether the symmetry condition (Eq. (19)
in Ref. [@sakuma-prb13]) on the unitary matrix
$\mathbf{U}^{(\mathbf{k})}$ is satisfied or not. See also Eq. (29) in
Ref. [@sakuma-prb13]. Used when `site_symmetry = .true`.

The default value is 1.0E-3.

#### `integer :: slwf_num`

The number of objective Wannier functions for selective localisation in
the selectively localised Wannier function (SLWF) method of
Ref. [@Marianetti]. These functions are obtained by minimising the
spread functional only with respect to the degrees of freedom of a
subset of `slwf_num` $<$ `num_wann` functions. At convergence, the
objective WFs will have a minimum cumulative spread, whereas the
remaining `num_wann` $-$ `slwf_num` functions are left unoptimised. The
initial guesses for the objective WFs are given by the first `slwf_num`
orbitals in the `projections` block. If `slwf_num = num_wann` no
selective minimisation is performed. In this case, `wannier90` will
simply generate a set of `num_wann` MLWFs.

The default is `num_wann`.

#### `logical :: slwf_constrain`

If `slwf_constrain=true`, then the centres of the objective Wannier
functions are constrained to either the centres of the first `slwf_num`
orbitals in the `projections` block or to new positions specified in the
`slwf_centres` block (see
Sec. [2.8.22](#sec:centre_constraints){reference-type="ref"
reference="sec:centre_constraints"}). In this case, a modified spread
functional, $\Omega_c$, with the addition of a constraint term, as
described in Ref. [@Marianetti].

The default is `false`

#### `real(kind=dp) :: slwf_lambda`

The value of the Lagrange multiplier $\lambda$ for the constraint term
in term added to modify the spread functional:
$\lambda \sum_{n=1}^{J'} \left(\overline{\mathbf{r}}_n - \mathbf{r}_{0n}\right)^2$,
where $J'$ is `slwf_num`, and $\overline{\mathbf{r}}_{n}$ and
$\mathbf{r}_{0n}$ are the centre and target centre, respectively, for
the $n^{\text{th}}$ objective WF.

The default is `0.0`.

#### Constraints on centres {#sec:centre_constraints}

If `slwf_constrain=true`, then by default the centres to which the
`slwf_num` objective Wannier function centres are constrained are given
by the first `slwf_num` rows of the `projections` block.

Optionally, the `slwf_centres` block may be used to define alternative
target centres for some or all of the `slwf_num` objective Wannier
functions.

The block below shows an example of how to set the constraints:

`begin slwf_centres`\
`   2  0.0   0.0  0.0`\
`   4  0.25  0.0  0.0`\
`end slwf_centres`

-   The first line sets the constraint for the centre of objective WF
    number 2 (as defined by the order of WFs in the `projections` block)
    to (0.0,0.0,0.0) in fractional co-ordinates.

-   The second line sets the constraint for the centre of objective WF
    number 4 (as defined by the order of WFs in the `projections` block)
    to (0.25,0.0,0.0) in fractional co-ordinates.

-   The target centres of all other objective Wannier functions remain
    as the centres given in the corresponding rows of the `projections`
    block.

### Post-Processing {#sec:post-p}

Capabilities:

-   Plot the WF

-   Plot the interpolated band structure

-   Plot the Fermi surface

-   Output the Hamiltonian in the WF basis

-   Transport calculation (quantum conductance and density of states)

#### `logical :: wannier_plot`

If $\verb#wannier_plot#=\verb#true#$, then the code will write out the
Wannier functions in a format specified by `wannier_plot_format`

The default value of this parameter is `false`.

#### `integer :: wannier_plot_list(:)`

A list of WF to plot. The WF numbered as per the `seedname.wout` file
after the minimisation of the spread.

The default behaviour is to plot all WF. For example, to plot WF 4, 5, 6
and 10:

`wannier_plot_list : 4-6, 10`

#### `integer :: wannier_plot_supercell`

The code generates the WFs on a grid corresponding to a
'super-unit-cell'. If `wannier_plot_supercell` is provided as a single
integer, then the size of the super-unit-cell is
`wannier_plot_supercell` times the size of the unit cell along all three
linear dimensions (the 'home' unit cell is kept approximately in the
middle); otherwise, if three integers are provided, the size of the
super-unit-cell is `wannier_plot_supercell(i)` times the size of the
unit cell along the $i-$th linear dimension.

The default value is 2.

#### `character(len=20) :: wannier_plot_format`

WF can be plotted in either XCrySDen (xsf) format or Gaussian cube
format. The valid options for this parameter are:

-   `xcrysden` (default)

-   `cube`

If `wannier_plot_format=xsf`: the code outputs the WF on the entire
super-unit-cell specified by `wannier_plot_supercell`.

If `wannier_plot_format=cube`: the code outputs the WF on a grid that is
smaller than the super-unit-cell specified by `wannier_plot_supercell`.
This grid is determined by `wannier_plot_mode`, `wannier_plot_radius`
and `wannier_plot_scale`, described in detail below.

The code is able to output Gaussian cube files for systems with
non-orthogonal lattice vectors. Many visualisation programs (including
XCrySDen), however, are only able to handle cube files for systems with
*orthogonal* lattice vectors. One visualisation program that is capable
of dealing with non-orthogonal lattice vectors is VESTA
(<http://jp-minerals.org/vesta/en/>).[^3]

#### `character(len=20) :: wannier_plot_mode`

Choose the mode in which to plot the WF, either as a molecule or as a
crystal.

The valid options for this parameter are:

-   `crystal` (default)

-   `molecule`

If `wannier_plot_format=cube`:

-   if `wannier_plot_mode = molecule`, then wherever the WF centre sits
    in the supercell, the origin of the cube is shifted (for the purpose
    of plotting only, ie, nothing is done to the U matrices etc) to
    coincide with the centre of mass of the atomic positions specified
    by the user in the `.win` input file. These atomic positions are
    also written to the cube file, so when it is visualised, the WF
    appears superimposed on the molecular structure.

-   if `wannier_plot_mode = crystal`, then the WF is not shifted, but
    instead the code searches for atoms that are within a radius of
    `wannier_plot_scale` $\times$ `wannier_plot_radius` of the WF centre
    and writes the coordinates of these atoms to the cube file. In this
    way, when the cube file is visualised, the WF appears superimposed
    on the nearest atoms to the WF centre.

-   `crystal` mode can be used for molecules, and `molecule` mode can be
    used for crystals.

#### `real(kind=dp) :: wannier_plot_radius`

If `wannier_plot_format=cube`, then ` wannier_plot_radius` is the radius
of the sphere that must fit inside the parallelepiped in which the WF is
plotted. `wannier_plot_radius` must be greater than 0. Units are Å.

The default value is 3.5.

#### `real(kind=dp) :: wannier_plot_scale`

If `wannier_plot_format=cube` and `wannier_plot_mode=crystal`, then the
code searches for atoms that are within a radius of `wannier_plot_scale`
$\times$ `wannier_plot_radius` of the WF centre and writes the
coordinates of these atoms to the cube file. In this way, when the cube
file is visualised, the WF appears superimposed on the nearest atoms to
the WF centre. `wannier_plot_scale` must be greater than 0. This
parameter is dimensionless.

The default value is 1.0.

#### `character(len=20) :: wannier_plot_spinor_mode`

If $\verb#spinors#=\verb#true#$ then this parameter controls the
quantity to plot. For a spinor WF with components $[\phi,\psi]$ the
quatity plotted is

-   `total` (default). $\sqrt{[|\phi|^2+|\psi|^2}$

-   `up`. $|\phi|\times sign(Re\{\phi\})$ if
    $\verb#wannier_plot_spinor_phase#=\verb#true#$, otherwise $|\phi|$

-   `down`. $|\psi|\times sign(Re\{\psi\})$ if
    $\verb#wannier_plot_spinor_phase#=\verb#true#$, otherwise $|\psi|$

Note: making a visual representation of a spinor WF is not as
straightforward as for a scalar WF. While a scalar WF is typically a
real valued function, a spinor WF is a complex, two component spinor.
`wannier90` is able to plot several different quantities derived from a
spinor WF which should give you a good idea of the nature of the WF.

#### `logical :: wannier_plot_spinor_phase`

If $\verb#wannier_plot_spinor_phase#=\verb#true#$ phase information will
be taken into account when plotting a spinor WF.

#### `logical :: bands_plot`

If $\verb#bands_plot#=\verb#true#$, then the code will calculate the
band structure, through Wannier interpolation, along the path in k-space
defined by `bands_kpath` using `bands_num_points` along the first
section of the path and write out an output file in a format specified
by `bands_plot_format`.

The default value is `false`.

#### kpoint_path

Defines the path in k-space along which to calculate the bandstructure.
Each line gives the start and end point (with labels) for a section of
the path. Values are in fractional coordinates with respect to the
primitive reciprocal lattice vectors.

`begin kpoint_path` $$\begin{array}{cccccccc}
G & 0.0 & 0.0 & 0.0 & L & 0.0 & 0.0 & 1.0 \\
L & 0.0 & 0.0 & 1.0 & N & 0.0 & 1.0 & 1.0 \\
\vdots
\end{array}$$ `end kpoint_path`

There is no default

#### `integer :: bands_num_points`

If $\verb#bands_plot#=\verb#true#$, then the number of points along the
first section of the bandstructure plot given by `kpoint_path`. Other
sections will have the same density of k-points.

The default value for `bands_num_points` is 100.

#### `character(len=20) :: bands_plot_format`

Format in which to plot the interpolated band structure. The valid
options for this parameter are:

-   `gnuplot` (default)

-   `xmgrace`

Note: it is possible to request output in both formats eg
$\verb#bands_format#=\verb#gnuplot xmgrace#$

#### `integer :: bands_plot_project(:)`

If present `wannier90` will compute the contribution of this set of WF
to the states at each point of the interpolated band structure. The WF
are numbered according to the seedname.wout file. The result is written
in the `seedname_band.dat` file, and a corresponding gnuplot script to
`seedname_band_proj.dat` .

For example, to project on to WFs 2, 6, 7, 8 and 12:

`bands_plot_project : 2, 6-8, 12`

#### `character(len=20) :: bands_plot_mode`

To interpolate the band structure along the k-point path, either use the
Slater-Koster interpolation scheme or truncate the Hamiltonian matrix in
the WF basis. Truncation criteria are provided by `hr_cutoff` and
`dist_cutoff`.

The valid options for this parameter are:

-   `s-k` (default)

-   `cut`

#### `integer :: bands_plot_dim`

Dimension of the system. If $\verb#bands_plot_dim#<\:$`<!-- -->`{=html}3
and $\verb#bands_plot_mode#=\verb#cut#$, lattice vector
$\mathbf{R}=N_1 \mathbf{A}_{1} + N_2 \mathbf{A}_{2} + N_3 \mathbf{A}_3$,
where $N_i=0$ if $\mathbf{A}_i$ is parallel to any of the confined
directions specified by `one_dim_axis`, are exclusively used in the band
structure interpolation.

The valid options for this parameter are:

-   3 (default)

-   2

-   1

#### `logical :: fermi_surface_plot`

If $\verb#fermi_surface_plot#=\verb#true#$, then the code will
calculate, through Wannier interpolation, the eigenvalues on a regular
grid with `fermi_surface_num_points` in each direction. The code will
write a file in bxsf format which can be read by XCrySDen in order to
plot the Fermi surface.

The default value is `false`.

#### `integer :: fermi_surface_num_points`

If $\verb#fermi_surface_plot#=\verb#true#$, then the number of divisions
in the regular k-point grid used to calculate the Fermi surface.

The default value for `fermi_surface_num_points` is 50.

#### `real(kind=dp) :: fermi_energy`

The Fermi energy in eV. This parameter is written into the bxsf file. If
`fermi_energy` is specified, ` fermi_energy_min`, `fermi_energy_max`,
and ` fermi_energy_step` should not be specified, and vice-versa.

The default value is 0.0

#### `real(kind=dp) :: fermi_energy_min`

Instead of specifyfing a single Fermi energy, it is possible to scan the
Fermi level over a range of values, and recompute certain quantities for
each $\varepsilon_F$.[^4] This is the minimum value in the range (in
eV).

There is no default value.

#### `real(kind=dp) :: fermi_energy_max`

The maximum value in the range of Fermi energies. Units are eV.

The default value is `fermi_energy_min`+1.0.

#### `real(kind=dp) :: fermi_energy_step`

Difference between consecutive values of the Fermi energy when scanning
from `fermi_energy_min` to ` fermi_energy_max`. Units are eV.

The default value is 0.01.

#### `character(len=20) :: fermi_surface_plot_format`

Format in which to plot the Fermi surface. The valid options for this
parameter are:

-   `xcrysden` (default)

#### `logical :: write_hr`

If $\verb#write_hr#=\verb#true#$, then the Hamiltonian matrix in the WF
basis will be written to a file `seedname_hr.dat`.

The default value is `false`.

#### `logical :: write_rmn`

If $\verb#write_rmn#=\verb#true#$, then the position operator in the WF
basis will be written to a file `seedname_r.dat`.

The default value is `false`.

#### `logical :: write_bvec`

If $\verb#write_bvec#=\verb#true#$, then the the matrix elements of
bvector and their weights will be written to a file `seedname.bvec`.

The default value is `false`.

#### `logical :: write_tb`

If $\verb#write_tb#=\verb#true#$, then the lattice vectors, together
with the Hamiltonian and position-operator matrices in the WF basis,
will be written to a file `seedname_tb.dat`, in units of Angstrom and
eV.

The default value is `false`.

#### `logical :: transport`

If $\verb#transport#=\verb#true#$, then the code will calculate quantum
conductance and density of states of a one-dimensional system. The
results will be written to files `seedname_qc.dat` and
`seedname_dos.dat`, respectively. Since both quantities are a function
of energy, they will be evaluated from `tran_win_min` to `tran_win_max`
with an interval of `tran_energy_step`.

The default value of this parameter is `false`.

#### `character(len=20) :: transport_mode`

If $\verb#transport_mode#=\verb#bulk#$, quantum conductance and density
of states are calculated for a perfectly-periodic one-dimensional
system. In this case, the transport part can either use the Hamiltonian
matrix in the WF basis generated by `wannier90` or a Hamiltonian matrix
provided by the external file `seedname_htB.dat`.

If $\verb#transport_mode#=\verb#lcr#$, quantum conductance and density
of states are calculated for a system where semi-infinite, left and
right leads are connected through a central conductor region. In this
case, the transport part will work independently from the
disentanglement and wannierise procedure. Details of the method is
described in Ref. [@nardelli-prb99].

If $\verb#tran_read_ht# = \verb#true#$ then the Hamiltonian matrices
must be provided by the five external files:
`seedname_htL.dat, seedname_htLC.dat, seedname_htC.dat, seedname_htCR.dat, seedname_htR.dat`.
If $\verb#tran_read_ht# = \verb#false#$ then the Hamiltonian matrices
are found automatically provided the supercell adheres to conditions
outlined in Section [7.3](#sec:2c2){reference-type="ref"
reference="sec:2c2"}.

The valid options for this parameter are:

-   `bulk` (default)

-   `lcr`

#### `real(kind=dp) :: tran_win_min`

The lower bound of the energy window for the transport calculation.
Units are eV.

The default value is -3.0.

#### `real(kind=dp) :: tran_win_max`

The upper bound of the energy window for the transport calculation.
Units are eV.

The default value is 3.0.

#### `real(kind=dp) :: tran_energy_step`

Sampling interval of the energy values from `tran_win_min` to
`tran_win_max`. Units are eV.

The default value is 0.01.

#### `real(kind=dp) :: fermi_energy`

The Fermi energy in eV. The energy axis of the quantum conductance and
density of states data will be shifted rigidly by this amount.

The default value is 0.0

#### `integer :: tran_num_bb`

Size of a bulk Hamiltonian matrix. This number is equal to the number of
WFs in one principal layer.

A one-dimensional system can be viewed as an array of principal layers
which are defined in a way that localized basis functions inside a
certain principal layer only interact with those in the nearest neighbor
principal layer. In `wannier90` a principal layer will be an integer
multiple of a unit cell, and the size is determined by `hr_cutoff`
and/or `dist_cutoff`. The criterion is rather arbitrary when WFs are
adopted as a localized basis set, and it is up to a user's choice.

The default value is 0.

#### `integer :: tran_num_ll`

Size of a left-lead Hamiltonian matrix. If
$\verb#transport_mode# = \verb#lcr#$ and
$\verb#tran_read_ht# = \verb#false#$ then `tran_num_ll` is the number of
Wannier functions in a principal layer.

The default value is 0.

#### `integer :: tran_num_rr`

Size of a right-lead Hamiltonian matrix.

The default value is 0.

#### `integer :: tran_num_cc`

Size of a conductor Hamiltonian matrix.

The default value is 0.

#### `integer :: tran_num_lc`

Number of columns in a left-lead_conductor Hamiltonian matrix. Number of
rows must be equal to `tran_num_ll`.

The default value is 0.

#### `integer :: tran_num_cr`

Number of rows in a conductor_right-lead Hamiltonian matrix. Number of
columns must be equal to `tran_num_rr`.

The default value is 0.

#### `integer :: tran_num_cell_ll`

Number of unit cells in one principal layer of left lead. Used if
$\verb#transport_mode# = \verb#lcr#$ and
$\verb#tran_read_ht# = \verb#false#$.

The default value is 0.

#### `integer :: tran_num_cell_rr`

Number of unit cells in one principal layer of right lead. Not used at
present.

The default value is 0.

#### `integer :: tran_num_bandc`

Half-bandwidth+1 of a band-diagonal conductor Hamiltonian matrix.

The Hamiltonian matrix of a central conductor part, which is read from
`seedname_htC.dat`, will be diagonally dominant when `tran_num_cc` is
very large. `tran_num_bandc` is used to construct a compact matrix which
contains the non-zero band-diagonal part of a full conductor Hamiltonian
matrix. Setting this parameter is only meaningful when `tran_num_bandc`
is greater than `tran_num_lc` and `tran_num_cr`.

The default value is 0.

#### `logical :: tran_write_ht`

If $\verb#tran_write_ht#=\verb#true#$, then the Hamiltonian matrix
formatted for the transport calculation will be written to a file
`seedname_htB.dat`.

The default value is `false`.

#### `logical :: tran_read_ht`

If $\verb#tran_write_ht#=\verb#true#$, then the Hamiltonian matrix
formatted for the transport calculation will be read from a set of files
described in the parameter `transport_mode`. Set
$\verb#tran_write_ht#=\verb#false#$ to perform automated lcr
calculations (see Section [7.3](#sec:2c2){reference-type="ref"
reference="sec:2c2"}).

The default value is `false`.

#### `logical :: tran_use_same_lead`

If $\verb#tran_use_same_lead#=\verb#true#$, then the left and the right
leads are the same. In this case, `seedname_htR.dat` is not required.

The default value is `true`.

#### `real(kind=dp) :: tran_group_threshold`

Used to group and sort Wannier functions according to the positions of
their centres. Wannier functions in a group are within
`tran_group_threshold` from one another in `x,y` and `z` directions.
Units are Å

The default is 0.15

#### `real(kind=dp) :: translation_centre_frac(3)`

Centre of the unit cell to which the final Wannier centres are
translated. Numbers are in fractional coordinates with respect to the
lattice vectors.

The default value is (0.0,0.0,0.0).

#### `logical :: use_ws_distance`

Improves the interpolation of the k-space Hamiltonian, by applying a
translation to each WF by a basis vector of the super-lattice that
minimises the distance between their centres. The translation is
dependent on both WF and on the unit cell vector to which they belong,
i.e., translate function $W_j({\bf r}-{\bf R})$ inside the Wigner-Seitz
cell centred on WF $W_i({\bf r})$.

For a longer explanation, see
Chapter [\[chap:interpolation\]](#chap:interpolation){reference-type="ref"
reference="chap:interpolation"}.

If `false` the code puts all the WF in the home cell, only possible
choice until wannier90 v2.0.1.

The default value is `true` (default changed since v.3.0). Introduced in
v2.1.

#### `real(kind=dp) :: ws_distance_tol`

Tolerance when determining whether two values
$\|\mathbf{d}_{ij\mathbf{R}} + \tilde{\mathbf{R}}_{nml} \|$ and
$\|\mathbf{d}_{ij\mathbf{R}} + \tilde{\mathbf{R}}_{n'm'l'} \|$ (as
defined in
chapter [\[chap:interpolation\]](#chap:interpolation){reference-type="ref"
reference="chap:interpolation"}) for the shortest distance between two
Wannier functions are equivalent. If the difference in distance (in
Angstrom) is less than `ws_distance_tol`, they are taken to be
equivalent.

The default value is $10^{-5}$.

#### `:: ws_search_size`

Maximum absolute value for the integers $n,m,l$ that identify the
super-lattice vectors $\tilde{\mathbf{R}}_{nml}$ (see
chapter [\[chap:interpolation\]](#chap:interpolation){reference-type="ref"
reference="chap:interpolation"}) when searching for points inside the
Wigner-Seitz cell. If `ws_search_size` is provided as a single integer,
then the number of repetitions of the Born-von Karman cell is the same
along all three linear dimensions; otherwise, if three integers are
provided, the number of repetitions along the $i-$th linear dimension is
`ws_search_size(i)`. The variable is used both in `hamiltonian.F90` and
in `ws_distance.F90`. In the latter case, its value is incremented by
one in order to account for WFs whose centre wanders away from the
original reference unit cell.\
The default value is generally sufficient, but might need to be
increased in case of elongated cells.

The default value is 2.

#### `logical :: write_u_matrices`

Write the $\mathbf{U}^{(\mathbf{k})}$ and
$\mathbf{U}^{\mathrm{dis}(\mathbf{k})}$ matrices obtained at the end of
wannierization to files `seedname_u.mat` and `seedname_u_dis.mat`,
respectively.

The default value is `false`.

#### `real(kind=dp) :: hr_cutoff`

The absolute value of the smallest matrix element of the Hamiltonian in
the WF basis. If $h_{mn}(\mathbf{R})>\:$` hr_cutoff`, then the matrix
element $h_{mn}(\mathbf{R})$ is retained and used in the band structure
interpolation (when $\verb#bands_plot_mode#=\verb#cut#$) or in the
transport calculation. Otherwise it is deemed to be insignificant and is
discarded. Units are eV.

The default value is 0.0.

#### `real(kind=dp) :: dist_cutoff`

The largest distance between two WFs for which the Hamiltonian matrix
element is retained and used in the band interpolation (when
$\verb#bands_plot_mode#=\verb#cut#$) or in the transport calculation.
Units are Å.

The default value is 1000.0.

#### `character(len=20) :: dist_cutoff_mode`

Dimension in which the distance between two WFs is calculated. The
vector connecting two WFs may be projected to a line (`one_dim`) or a
plane (`two_dim`). The size of the projected vector is calculated, and
`dist_cutoff` is applied. When `one_dim` or `two_dim` is used,
`one_dim_axis` must be given to specify extended or confined direction.

The valid options for this parameter are:

-   `three_dim` (default)

-   `two_dim`

-   `one_dim`

#### `character(len=20) :: one_dim_axis`

Extended direction for a one-dimensional system or confined direction
for a two-dimensional system. This direction must be parallel to one of
the Cartesian axes.

The valid options for this parameter are:

-   `x`

-   `y`

-   `z`

No default.

## Projections {#ch:proj}

### Specification of projections in `seedname.win` {#sec:proj}

Here we describe the projection functions used to construct the initial
guess $A_{mn}^{(\mathbf{k})}$ for the unitary transformations.

Each projection is associated with a site and an angular momentum state
defining the projection function. Optionally, one may define, for each
projection, the spatial orientation, the radial part, the diffusivity,
and the volume over which real-space overlaps $A_{mn}$ are calculated.

The code is able to

1.  project onto s,p,d and f angular momentum states, plus the hybrids
    sp, sp$^2$, sp$^3$, sp$^3$d, sp$^3$d$^2$.

2.  control the radial part of the projection functions to allow higher
    angular momentum states, e.g., both 3s and 4s in silicon.

The atomic orbitals of the hydrogen atom provide a good basis to use for
constructing the projection functions: analytical mathematical forms
exist in terms of the good quantum numbers $n$, $l$ and $m$; hybrid
orbitals (sp, sp$^{2}$, sp$^{3}$, sp$^{3}$d etc.) can be constructed by
simple linear combination $|\phi\rangle =
\sum_{nlm} C_{nlm}|nlm\rangle$ for some coefficients $C_{nlm}$.

The angular functions that use as a basis for the projections are not
the canonical spherical harmonics $Y_{lm}$ of the hydrogenic Schrödinger
equation but rather the *real* (in the sense of non-imaginary) states
$\Theta_{lm_{\mathrm{r}}}$, obtained by a unitary transformation. For
example, the canonical eigenstates associated with $l=1$, $m=\{-1,0,1\}$
are not the real p$_{x}$, p$_{y}$ and p$_{z}$ that we want. See
Section [3.4](#sec:orbital-defs){reference-type="ref"
reference="sec:orbital-defs"} for our mathematical conventions regarding
projection orbitals for different $n$, $l$ and $m_{\mathrm{r}}$.

We use the following format to specify projections in `<seedname>.win`:

`Begin Projections`\
`[units]`\
`site:ang_mtm:zaxis:xaxis:radial:zona`\
`    `⋮\
`End Projections`

Notes:

`units`:\
Optional. Either `Ang` or `Bohr` to specify whether the projection
centres specified in this block (if given in Cartesian co-ordinates) are
in units of Angstrom or Bohr, respectively. The default value is `Ang`.

`site`:\
`C`, `Al`, etc. applies to all atoms of that type\
`f=0,0.50,0` -- centre on (0.0,0.5,0.0) in **f**ractional coordinates
(crystallographic units) relative to the direct lattice vectors\
`c=0.0,0.805,0.0` -- centre on (0.0,0.805,0.0) in **C**artesian
coordinates in units specified by the optional string `units` in the
first line of the projections block (see above).

`ang_mtm`:\
Angular momentum states may be specified by `l` and `mr`, or by the
appropriate character string. See
Tables [3.1](#tab:angular){reference-type="ref" reference="tab:angular"}
and [3.2](#tab:hybrids){reference-type="ref" reference="tab:hybrids"}.
Examples:\
`l=2,mr=1 ` or ` dz2` -- a single projection with $l=2$,
$m_{\textrm{r}}=1$ (i.e., d$_{z^{2}}$)\
`l=2,mr=1,4 ` or ` dz2,dx2-y2` -- two functions: d$_{z^{2}}$ and
d$_{xz}$\
`l=-3 ` or ` sp3` -- four sp$^{3}$ hybrids\
Specific hybrid orbitals may be specified as follows:\
`l=-3,mr=1,3 ` or ` sp3-1,sp3-3` -- two specific sp$^{3}$ hybrids\
Multiple states may be specified by separating with '`;`', e.g.,\
`sp3;l=0 ` or ` l=-3;l=0` -- four sp$^{3}$ hybrids and one s orbital

`zaxis` (optional):\
`z=1,1,1` -- set the $z$-axis to be in the (1,1,1) direction. Default is
`z=0,0,1`

`xaxis` (optional):\
`x=1,1,1` -- set the $x$-axis to be in the (1,1,1) direction. Default is
`x=1,0,0`

`radial` (optional):\
`r=2` -- use a radial function with one node (ie second highest
pseudostate with that angular momentum). Default is `r=1`. Radial
functions associated with different values of `r` should be orthogonal
to each other.

`zona` (optional):\
`zona=2.0` -- the value of $\frac{Z}{a}$ for the radial part of the
atomic orbital (controls the diffusivity of the radial function). Units
always in reciprocal Angstrom. Default is `zona=1.0`.

**Examples**

1\. CuO, s,p and d on all Cu; sp$^3$ hybrids on O:

`Cu:l=0;l=1;l=2 `

`O:l=-3 ` or ` O:sp3`

2\. A single projection onto a p$_z$ orbital orientated in the (1,1,1)
direction:

`c=0,0,0:l=1,mr=1:z=1,1,1 ` or ` c=0,0,0:pz:z=1,1,1`

3\. Project onto s, p and d (with no radial nodes), and s and p (with
one radial node) in silicon:

`Si:l=0;l=1;l=2`

`Si:l=0;l=1:r=2`

### Spinor Projections

When `spinors=.true.` it is possible to select a set of localised
functions to project onto 'up' states and a set to project onto 'down'
states where, for complete flexibility, it is also possible to set the
local spin quantisation axis.

Note, however, that this feature requires a recent version of the
interface between the ab-initio code and Wannier90 (i.e., written after
the release of the 2.0 version, in October 2013) supporting spinor
projections.

`Begin Projections`\
`[units]`\
`site:ang_mtm:zaxis:xaxis:radial:zona(spin)[quant_dir]`\
`    `⋮\
`End Projections`

`spin` (optional):\
Choose projection onto 'up' or 'down' states\
`u` -- project onto 'up' states.\
`d` -- project onto 'down' states.\
Default is `u,d`

`quant_dir` (optional):\
`1,0,0` -- set the spin quantisation axis to be in the (1,0,0)
direction. Default is `0,0,1`

**Examples**

-   18 projections on an iron site

    `Fe:sp3d2;dxy;dxx;dyz`

-   same as above

    `Fe:sp3d2;dxy;dxx;dyz(u,d)`

-   same as above

    `Fe:sp3d2;dxy;dxz;dyz(u,d)[0,0,1]`

-   same as above but quantisation axis is now x

    `Fe:sp3d2;dxy;dxz;dyz(u,d)[1,0,0]`

-   now only 9 projections onto up states

    `Fe:sp3d2;dxy;dxz;dyz(u)`

-   9 projections onto up-states and 3 on down

    `Fe:sp3d2;dxy;dxz;dyz(u) `\
    `Fe:dxy;dxz;dyz(d)`

-   projections onto alternate spin states for two lattice sites (Cr1,
    Cr2)

    `Cr1:d(u)`\
    `Cr2:d(d)`

### Short-Cuts

#### Random projections

It is possible to specify the projections, for example, as follows:

`Begin Projections`\
`random`\
`C:sp3`\
`End Projections`

in which case `wannier90` uses four sp$^3$ orbitals centred on each C
atom and then chooses the appropriate number of randomly-centred s-type
Gaussian functions for the remaining projection functions. If the block
only consists of the string `random` and no specific projection centres
are given, then all of the projection centres are chosen randomly.

#### Bloch phases

Setting `use_bloch_phases = true` in the input file absolves the user of
the need to specify explicit projections. In this case, the Bloch
wave-functions are used as the projection orbitals, namely
$A_{mn}^{(\mathbf{k})} =
\langle\psi_{m\mathbf{k}}|\psi_{n\mathbf{k}}\rangle = \delta_{mn}$.

### Orbital Definitions {#sec:orbital-defs}

The angular functions $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$
associated with particular values of $l$ and $m_{\mathrm{r}}$ are given
in Tables [3.1](#tab:angular){reference-type="ref"
reference="tab:angular"} and [3.2](#tab:hybrids){reference-type="ref"
reference="tab:hybrids"}.

The radial functions $R_{\mathrm{r}}(r)$ associated with different
values of $r$ should be orthogonal. One choice would be to take the set
of solutions to the radial part of the hydrogenic Schrödinger equation
for $l=0$, i.e., the radial parts of the 1s, 2s, 3s... orbitals, which
are given in Table [3.3](#tab:radial){reference-type="ref"
reference="tab:radial"}.

::: center
::: {#tab:angular}
  ----- ------------------ -------------- ---------------------------------------------------------------------------------------------
                                          
   $l$   $m_{\mathrm{r}}$       Name                               $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$
                                          
                                          
    0           1               `s`                                          $\frac{1}{\sqrt{4\pi}}$
                                          
                                          
    1           1               `pz`                                    $\sqrt{\frac{3}{4\pi}}\cos\theta$
                                          
    1           2               `px`                              $\sqrt{\frac{3}{4\pi}}\sin\theta\cos\varphi$
                                          
    1           3               `py`                              $\sqrt{\frac{3}{4\pi}}\sin\theta\sin\varphi$
                                          
                                          
    2           1              `dz2`                              $\sqrt{\frac{5}{16\pi}}(3\cos^{2}\theta -1)$
                                          
    2           2              `dxz`                         $\sqrt{\frac{15}{4\pi}}\sin\theta\cos\theta\cos\varphi$
                                          
    2           3              `dyz`                         $\sqrt{\frac{15}{4\pi}}\sin\theta\cos\theta\sin\varphi$
                                          
    2           4             `dx2-y2`                         $\sqrt{\frac{15}{16\pi}}\sin^{2}\theta\cos2\varphi$
                                          
    2           5              `dxy`                           $\sqrt{\frac{15}{16\pi}}\sin^{2}\theta\sin2\varphi$
                                          
                                          
    3           1              `fz3`                       $\frac{\sqrt{7}}{4\sqrt{\pi}}(5\cos^{3}\theta-3\cos\theta)$
                                          
    3           2              `fxz2`               $\frac{\sqrt{21}}{4\sqrt{2\pi}}(5\cos^{2}\theta-1)\sin\theta\cos\varphi$
                                          
    3           3              `fyz2`               $\frac{\sqrt{21}}{4\sqrt{2\pi}}(5\cos^{2}\theta-1)\sin\theta\sin\varphi$
                                          
    3           4           `fz(x2-y2)`               $\frac{\sqrt{105}}{4\sqrt{\pi}}\sin^{2}\theta\cos\theta\cos2\varphi$
                                          
    3           5              `fxyz`                 $\frac{\sqrt{105}}{4\sqrt{\pi}}\sin^{2}\theta\cos\theta\sin2\varphi$
                                          
    3           6           `fx(x2-3y2)`   $\frac{\sqrt{35}}{4\sqrt{2\pi}}\sin^{3}\theta(\cos^{2}\varphi-3\sin^{2}\varphi)\cos\varphi$
                                          
    3           7           `fy(3x2-y2)`   $\frac{\sqrt{35}}{4\sqrt{2\pi}}\sin^{3}\theta(3\cos^{2}\varphi-\sin^{2}\varphi)\sin\varphi$
                                          
  ----- ------------------ -------------- ---------------------------------------------------------------------------------------------

  : Angular functions $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$
  associated with particular values of $l$ and $m_{\mathrm{r}}$ for
  $l\ge0$.
:::
:::

::: center
::: {#tab:hybrids}
  ----------------------- ------------------ ----------- -----------------------------------------------------------------------------
                                                         
            $l$            $m_{\mathrm{r}}$     Name                      $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$
                                                         
                                                         
   $-$`<!-- -->`{=html}1          1            `sp-1`                  $\frac{1}{\sqrt{2}}$`s` $+\frac{1}{\sqrt{2}}$`px`
                                                         
   $-$`<!-- -->`{=html}1          2            `sp-2`                  $\frac{1}{\sqrt{2}}$`s` $-\frac{1}{\sqrt{2}}$`px`
                                                         
                                                         
   $-$`<!-- -->`{=html}2          1            `sp2-1`    $\frac{1}{\sqrt{3}}$`s` $-\frac{1}{\sqrt{6}}$`px` $+\frac{1}{\sqrt{2}}$`py`
                                                         
   $-$`<!-- -->`{=html}2          2            `sp2-2`    $\frac{1}{\sqrt{3}}$`s` $-\frac{1}{\sqrt{6}}$`px` $-\frac{1}{\sqrt{2}}$`py`
                                                         
   $-$`<!-- -->`{=html}2          3            `sp2-3`                 $\frac{1}{\sqrt{3}}$`s` $+\frac{2}{\sqrt{6}}$`px`
                                                         
                                                         
   $-$`<!-- -->`{=html}3          1            `sp3-1`                   $\frac{1}{2}$(`s` $+$ `px` $+$ `py` $+$ `pz`)
                                                         
   $-$`<!-- -->`{=html}3          2            `sp3-2`                   $\frac{1}{2}$(`s` $+$ `px` $-$ `py` $-$ `pz`)
                                                         
   $-$`<!-- -->`{=html}3          3            `sp3-3`                   $\frac{1}{2}$(`s` $-$ `px` $+$ `py` $-$ `pz`)
                                                         
   $-$`<!-- -->`{=html}3          4            `sp3-4`                   $\frac{1}{2}$(`s` $-$ `px` $-$ `py` $+$ `pz`)
                                                         
                                                         
   $-$`<!-- -->`{=html}4          1           `sp3d-1`    $\frac{1}{\sqrt{3}}$`s` $-\frac{1}{\sqrt{6}}$`px` $+\frac{1}{\sqrt{2}}$`py`
                                                         
   $-$`<!-- -->`{=html}4          2           `sp3d-2`    $\frac{1}{\sqrt{3}}$`s` $-\frac{1}{\sqrt{6}}$`px` $-\frac{1}{\sqrt{2}}$`py`
                                                         
   $-$`<!-- -->`{=html}4          3           `sp3d-3`                 $\frac{1}{\sqrt{3}}$`s` $+\frac{2}{\sqrt{6}}$`px`
                                                         
   $-$`<!-- -->`{=html}4          4           `sp3d-4`                $\frac{1}{\sqrt{2}}$`pz` $+\frac{1}{\sqrt{2}}$`dz2`
                                                         
   $-$`<!-- -->`{=html}4          5           `sp3d-5`               $-\frac{1}{\sqrt{2}}$`pz` $+\frac{1}{\sqrt{2}}$`dz2`
                                                         
                                                         
   $-$`<!-- -->`{=html}5          1           `sp3d2-1`             $\frac{1}{\sqrt{6}}\verb#s#-\frac{1}{\sqrt{2}}\verb#px#
                                                                   -\frac{1}{\sqrt{12}}\verb#dz2#+\frac{1}{2}\verb#dx2-y2#$
                                                         
   $-$`<!-- -->`{=html}5          2           `sp3d2-2`             $\frac{1}{\sqrt{6}}\verb#s#+\frac{1}{\sqrt{2}}\verb#px#
                                                                   -\frac{1}{\sqrt{12}}\verb#dz2#+\frac{1}{2}\verb#dx2-y2#$
                                                         
   $-$`<!-- -->`{=html}5          3           `sp3d2-3`             $\frac{1}{\sqrt{6}}\verb#s#-\frac{1}{\sqrt{2}}\verb#py#
                                                                   -\frac{1}{\sqrt{12}}\verb#dz2#-\frac{1}{2}\verb#dx2-y2#$
                                                         
   $-$`<!-- -->`{=html}5          4           `sp3d2-4`             $\frac{1}{\sqrt{6}}\verb#s#+\frac{1}{\sqrt{2}}\verb#py#
                                                                   -\frac{1}{\sqrt{12}}\verb#dz2#-\frac{1}{2}\verb#dx2-y2#$
                                                         
   $-$`<!-- -->`{=html}5          5           `sp3d2-5`             $\frac{1}{\sqrt{6}}\verb#s#-\frac{1}{\sqrt{2}}\verb#pz#
                                                                                +\frac{1}{\sqrt{3}}\verb#dz2#$
                                                         
   $-$`<!-- -->`{=html}5          6           `sp3d2-6`             $\frac{1}{\sqrt{6}}\verb#s#+\frac{1}{\sqrt{2}}\verb#pz#
                                                                                +\frac{1}{\sqrt{3}}\verb#dz2#$
                                                         
  ----------------------- ------------------ ----------- -----------------------------------------------------------------------------

  : Angular functions $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$
  associated with particular values of $l$ and $m_{\mathrm{r}}$ for
  $l<0$, in terms of the orbitals defined in
  Table [3.1](#tab:angular){reference-type="ref"
  reference="tab:angular"}.
:::
:::

::: center
::: {#tab:radial}
  ---------- --------------------------------------------
             
     $r$                 $R_{\mathrm{r}}(r)$
             
             
      1            $2 \alpha^{3/2}\exp(-\alpha r)$
             
             
      2       $\frac{1}{2\sqrt{2}}\alpha^{3/2}(2-\alpha
                         r)\exp(-\alpha r/2)$
             
             
      3       $\sqrt{\frac{4}{27}}\alpha^{3/2}(1-2\alpha
              r/3+2\alpha^{2}r^{2}/27)\exp(-\alpha r/3)$
             
  ---------- --------------------------------------------

  :  One possible choice for the radial functions $R_{\mathrm{r}}(r)$
  associated with different values of $r$: the set of solutions to the
  radial part of the hydrogenic Schrödinger equation for $l=0$, i.e.,
  the radial parts of the 1s, 2s, 3s... orbitals, where
  $\alpha=Z/a={\tt zona}$.
:::
:::

### Projections via the SCDM-**k** method in pw2wannier90

For many systems, such as aperiodic systems, crystals with defects, or
novel materials with complex band structure, it may be extremely hard to
identify *a-priori* a good initial guess for the projection functions
used to generate the $A_{mn}^{(\mathbf{k})}$ matrices. In these cases,
one can use a different approach, known as the SCDM-**k**
method[@LinLin-ArXiv2017], based on a QR factorization with column
pivoting (QRCP) of the density matrix from the self-consistent field
calculation, which allows one to avoid the tedious step of specifying a
projection block altogether, hence to avoid . This method is robust in
generating well localised function with the correct spatial orientations
and in general in finding the global minimum of the spread functional
$\Omega$. Any electronic-structure code should in principle be able to
implement the SCDM-**k** method within their interface with Wannier90,
however at the moment (develop branch on the GitHub repository July
2019) only the Quantum ESPRESSO package has this capability implemented
in the `pw2wannier90` interface program. Moreover, the `pw2wannier90`
interface program supports also the SCDM-**k** method for
spin-noncollinear systems. The SCDM-**k** can operate in two modes:

1.  In isolation, i.e., without performing a subsequent Wannier90
    optimisation (not recommended). This can be achieved by setting
    `num_iter=0` and `dis_num_iter=0` in the `<seedname>.win` input
    file. The rationale behind this is that in general the projection
    functions obtained with the SCDM-**k** are already well localised
    with the correct spatial orientations. However, the spreads of the
    resulting functions are usually larger than the MLWFs ones.

2.  In combination with the Marzari-Vanderbilt (recommended option). In
    this case, the SCDM-**k** is only used to generate the initial
    $A_{mn}^{(\mathbf{k})}$ matrices as a replacement scheme for the
    projection block.

The following keywords need to be specified in the `pw2wannier90.x`
input file `<seedname>.pw2wan`: `scdm_proj` `scdm_entanglement`
`scdm_mu` `scdm_sigma`

### Projections via pseudo-atomic orbitals in pw2wannier90 {#sec:proj_pdwf}

When generating pseudopotentials, often the atomic wavefunctions of
isolated atom are pseudized and bundled together with the
pseudopotential files. These orbitals are often used for computing the
projectabilities, for instance, measuring orbital contributions to band
structures. Instead of manually specifying the initial projections in
the `projections` block, one can use these pseudo-atomic orbitals to
automate the initial projection process.

Currently (July 2023), this functionality is implemented in the
[quantum-espresso]{.smallcaps} interface, but in principle it can be
done in any other interface as well. In the following, we will use the
[quantum-espresso]{.smallcaps} interface as an example to illustrate the
whole procedure.

To activate pseudo-atomic orbital projection, one needs to set
`auto_projections = .true.` in the `win` file, and remove the
`projections` block.

Then in the `pw2wannier90` input file, one needs to add an additional
tag `atom_proj = .true.`. This will ask `pw2wannier90` to read the
pseudo-atomic orbitals from the pseudopotential files, and use them to
compute the `amn` file.

Some times, one may want to exclude semi-core states from
Wannierisation, for such cases, one can inspect the stdout of
`pw2wannier90`, which will print the orbitals used for computing `amn`,
e.g.,

` `

> -------------------------------------\
> \*\*\* Compute A with atomic projectors\
> -------------------------------------\
> Use atomic projectors from UPF\
> \
> (read from pseudopotential files):\
> state \# 1: atom 1 (C ), wfc 1 (l=0 m= 1)\
> state \# 2: atom 1 (C ), wfc 2 (l=1 m= 1)\
> state \# 3: atom 1 (C ), wfc 2 (l=1 m= 2)\
> state \# 4: atom 1 (C ), wfc 2 (l=1 m= 3)\
> state \# 5: atom 2 (C ), wfc 1 (l=0 m= 1)\
> state \# 6: atom 2 (C ), wfc 2 (l=1 m= 1)\
> state \# 7: atom 2 (C ), wfc 2 (l=1 m= 2)\
> state \# 8: atom 2 (C ), wfc 2 (l=1 m= 3)\

Here it shows that there are two carbon atoms, each with one $s$ and
three $p$ orbitals. If one wants to exclude specific orbital(s), there
is an additional input `atom_proj_exclude`, which accept a list of
integers, e.g.,

`atom_proj_exclude = 1 5`

which will exclude the two $s$ orbitals from computing `amn`.

##### Advanced usage {#advanced-usage .unnumbered}

If the pseudopotential orbitals are not enough, one could also generate
a custom set of orbitals, and ask `pw2wannier90` to use them for
computing `amn`. This can be done by setting

`atom_proj_dir = ’./ext_proj’`

where the directory `ext_proj` contains the orbitals for all the atomic
species used in the calculation. For example, for a silicon calculation,
the directory `ext_proj` should contain a file named `Si.dat`. The
format of the file is:

1.  The first line contains two integers: the number of radial grid
    points ($n_g$) and the number of projectors ($n_p$), e.g.,

        1141 2

    which means the radial grid has $n_g = 1141$ points, and there are
    $n_p = 3$ projectors.

2.  The second line contains $n_p$ integers specifying the angular
    momentums of all the projectors, e.g.,

        0 1

    standing for the two projectors having $s$ and $p$ characters,
    respectively.

3.  The rest of the file contains $n_g$ rows of the radial wavefunctions
    of the projectors. There are $2+n_p$ columns: the first column is
    the $x$-grid, the second column is the $r$-grid in Bohr unit, and
    they are related by $r = \exp(x)$. The rest are $n_p$ columns of the
    radial wavefunctions of the projectors,

        -9.639057329615259 0.000065134426111 3.32211124436945e-05 1.86840239681223e-09
        -9.626557329615258 0.000065953716334 3.363898259696903e-05 1.915701228607072e-09
        -9.614057329615258 0.000066783311958 3.406210890972733e-05 1.964197436025957e-09
        ...

    Inside `pw2wannier90.x`, the radial wavefunction will be read and
    multiplied by spherical harmonics to form the actual projectors.

    For a practical example of extracting pseudo-atomic orbitals from
    UPF file and writing to a `pw2wannier90`-recognizable `.dat` file,
    see the script `utility/write_pdwf_projectors.py`.

    For an actual example of a `Si.dat` file for silicon, see the file
    `examples/example35/ext_proj/Si.dat`.

## Code Overview

`wannier90` can operate in two modes:

1.  *Post-processing mode:* read in the overlaps and projections from
    file as computed inside a first-principles code. We expect this to
    be the most common route to using `wannier90`, and is described in
    Ch. [5](#ch:wann-pp){reference-type="ref" reference="ch:wann-pp"};

2.  *Library mode:* as a set of library routines to be called from
    within a first-principles code that passes the overlaps and
    projections to the `wannier90` library routines and in return gets
    the unitary transformation corresponding to MLWF. This route should
    be used if the MLWF are needed within the first-principles code, for
    example in post-LDA methods such as LDA+U or SIC, and is described
    in Ch. [6](#ch:wann-lib){reference-type="ref"
    reference="ch:wann-lib"}.

<figure id="structure">
<div class="center">
<embed src="overview.pdf" style="width:5in" />
</div>
<figcaption>Schematic overview of the module structure of
<code>wannier90</code>. Modules may only use data and subroutines from
lower modules.</figcaption>
</figure>

## `wannier90` as a post-processing tool {#ch:wann-pp}

This is a description of how to use `wannier90` as a post-processing
tool.

The code must be run twice. On the first pass either the logical keyword
`postproc_setup` must be set to `.true.` in the input file
`seedname.win` or the code must be run with the command line option
`-pp`. Running the code then generates the file `seedname.nnkp` which
provides the information required to construct the
$M_{mn}^{(\mathbf{k,b})}$ overlaps (Ref. [@marzari-prb97], Eq. (25)) and
$A_{mn}^{(\mathbf{k})}$ (Ref. [@marzari-prb97], Eq. (62);
Ref. [@souza-prb01], Eq. (22)).

Once the overlaps and projection have been computed and written to files
`seedname.mmn` and `seedname.amn`, respectively, set `postproc_setup` to
`.false.` and run the code. Output is written to the file
`seedname.wout`.

### `seedname.nnkp` file

OUTPUT, if $\verb#postproc_setup#=\verb#.true.#$

The file `seedname.nnkp` provides the information needed to determine
the required overlap elements $M_{mn}^{(\mathbf{k,b})}$ and projections
$A_{mn}^{(\mathbf{k})}$. It is written automatically when the code is
invoked with the `-pp` command-line option (or when
`postproc_setup=.true.` in `seedname.win`. There should be no need for
the user to edit this file.

Much of the information in `seedname.nnkp` is arranged in blocks
delimited by the strings `begin block_name` ... `end block_name`, as
described below.

#### Keywords

The first line of the file is a user comment, e.g., the date and time:

`File written on 12Feb2006 at 15:13:12`

The only logical keyword is `calc_only_A`, eg,

`calc_only_A  :  F`

#### `Real_lattice` block

    begin real_lattice
     2.250000   0.000000   0.000000
     0.000000   2.250000   0.000000
     0.000000   0.000000   2.250000
    end real_lattice

The real lattice vectors in units of Angstrom.

#### `Recip_lattice` block

    begin recip_lattice
     2.792527   0.000000   0.000000
     0.000000   2.792527   0.000000
     0.000000   0.000000   2.792527
    end recip_lattice

The reciprocal lattice vectors in units of inverse Angstrom.

#### `Kpoints` block

    begin kpoints
      8
      0.00000   0.00000   0.00000
      0.00000   0.50000   0.00000
      .
      .
      .
      0.50000   0.50000   0.50000
    end kpoints

The first line in the block is the total number of k-points `num_kpts`.
The subsequent `num_kpts` lines specify the k-points in crystallographic
co-ordinates relative to the reciprocal lattice vectors.

#### `Projections` block

    begin projections
       n_proj
       centre   l  mr  r   
         z-axis   x-axis   zona
       centre   l  mr  r   
         z-axis   x-axis   zona
       .
       .
    end projections

Notes:

`n_proj`: integer; the number of projection centres, equal to the number
of MLWF `num_wann`.

`centre`: three real numbers; projection function centre in
crystallographic co-ordinates relative to the direct lattice vectors.

`l  mr  r`: three integers; $l$ and $m_\mathrm{r}$ specify the angular
part $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$, and $\mathrm{r}$
specifies the radial part $R_{\mathrm{r}}(r)$ of the projection function
(see Tables [3.1](#tab:angular){reference-type="ref"
reference="tab:angular"}, [3.2](#tab:hybrids){reference-type="ref"
reference="tab:hybrids"} and [3.3](#tab:radial){reference-type="ref"
reference="tab:radial"}).

`z-axis`: three real numbers; default is `0.0 0.0 1.0`; defines the axis
from which the polar angle $\theta$ in spherical polar coordinates is
measured.

`x-axis`: three real numbers; must be orthogonal to `z-axis`; default is
`1.0 0.0 0.0` or a vector perpendicular to `z-axis` if `z-axis` is
given; defines the axis from with the azimuthal angle $\varphi$ in
spherical polar coordinates is measured.

`zona`: real number; the value of $\frac{Z}{a}$ associated with the
radial part of the atomic orbital. Units are in reciprocal Angstrom.

#### `spinor_projections` block

    begin spinor_projections
       n_proj
       centre   l  mr  r   
        z-axis   x-axis   zona
         spin spn_quant
       centre   l  mr  r   
        z-axis   x-axis   zona
         spin spn_quant
       .
       .
    end spinor_projections

Notes: Only one of projections and spinor_projections should be defined.
Variables are the same as the projections block with the addition of
`spin` and `spn_quant`.

`spin`: integer. '1' or '-1' to denote projection onto up or down
states.

`spn_quant`: three real numbers. Defines the spin quantisation axis in
Cartesian coordinates.

#### `nnkpts` block

    begin nnkpts
      10
      1   2   0  0  0
      .
      .
    end nnkpts

First line: `nntot`, the number of nearest neighbours belonging to each
k-point of the Monkhorst-Pack mesh

Subsequent lines: `nntot`$\times$`num_kpts` lines, ie, `nntot` lines of
data for each k-point of the mesh.

Each line of consists of 5 integers. The first is the k-point number
`nkp`. The second to the fifth specify it's nearest neighbours
$\mathbf{k+b}$: the second integer points to the k-point that is the
periodic image of the $\mathbf{k+b}$ that we want; the last three
integers give the G-vector, in reciprocal lattice units, that brings the
k-point specified by the second integer (which is in the first BZ) to
the actual $\mathbf{k+b}$ that we need.

#### `exclude_bands` block

    begin exclude_bands 
      8 
      1 
      2 
      .
      .
    end exclude_bands

To exclude bands (independent of k-point) from the calculation of the
overlap and projection matrices, for example to ignore shallow-core
states. The first line is the number of states to exclude, the following
lines give the states for be excluded.

#### []{#sec:auto-projections-block label="sec:auto-projections-block"}`auto_projections` block

    begin auto_projections
       8
       0
    end auto_projections

This block is only printed if `auto_projections=true` in the input. The
choice of an additional block has been made in order to maintain
back-compatibility with codes that interface with `wannier90`, e.g.
`pw2wannier90`. The first entry in the block (in the example above, `8`)
is the total number of target projections and it is equal to the number
of sought Wannier functions.

The second entry is a reserved flag with the value of zero. The
implementations of the interface codes MUST check for this value to be
zero and stop otherwise. In the future, one possible extension that we
plan is to combine the automatic generation of initial projections with
the selection of projections via a projections block. This will allow
the user to specify only a subset of initial projections in the
projections block and leave the interface code to automatically generate
the remaining ones. In that case the constraint on the second entry will
be lifted, so that it can take on the meaning of the number of
projections that need to be generated automatically.

The selected columns of the density matrix (SCDM)
method [@LinLin-ArXiv2017] is one way of generating the initial
$A_{mn}^{(\mathbf{k})}$ in an automatic way. This has been implemented
in the `pw2wannier90` interface code (you need v6.3 with the files
provided in the `pwscf` folder of Wannier90, or v6.4), see for instance
Example 27 in the `wannier90` tutorial that shows how to use it.

Moreover, also the automatic generation of initial projections with
spinor WFs is implemented in the `pw2wannier90` interface. See Example
31 in the `wannier90` tutorial that shows how to use it.

Another automatic projection method is projectability-disentangled
Wannier function (PDWF) [@Qiao2023-pdwf], which uses pseudo-atomic
orbitals inside pseudopotentials as initial guesses. See Example 34 and
35.

#### An example of projections {#sec:proj_example}

As a concrete example: one wishes to have a set of four sp$^3$
projection orbitals on, say, a carbon atom at (0.5,0.5,0.5) in
fractional co-ordinates relative to the direct lattice vectors. In this
case `seedname.win` will contain the following lines:

    begin projections
     C:l=-1
    end projections

and `seedname.nnkp`, generated on the first pass of `wannier90` (with
`postproc_setup=T`), will contain:

    begin projections
       4
       0.50000    0.50000    0.50000    -1  1  1
         0.000  0.000  1.000   1.000  0.000  0.000   2.00 
       0.50000    0.50000    0.50000    -1  2  1
         0.000  0.000  1.000   1.000  0.000  0.000   2.00 
       0.50000    0.50000    0.50000    -1  3  1
         0.000  0.000  1.000   1.000  0.000  0.000   2.00 
       0.50000    0.50000    0.50000    -1  4  1
         0.000  0.000  1.000   1.000  0.000  0.000   2.00 
    end projections

where the first line tells us that in total four projections are
specified, and the subsquent lines provide the projection centre, the
angular and radial parts of the orbital (see
Section [3.4](#sec:orbital-defs){reference-type="ref"
reference="sec:orbital-defs"} for definitions), the $z$ and $x$ axes,
and the diffusivity and cut-off radius for the projection orbital.

[pwscf]{.smallcaps}, or any other *ab initio* electronic structure code,
then reads `seedname.nnkp` file, calculates the projections and writes
them to `seedname.amn`.

### `seedname.mmn` file

INPUT.

The file `seedname.mmn` contains the overlaps $M_{mn}^{(\mathbf{k,b})}$.

First line: a user comment, e.g., the date and time

Second line: 3 integers: `num_bands`, `num_kpts`, `nntot`

Then: $\verb#num_kpts#\times\verb#nntot#$ blocks of data:

First line of each block: 5 integers. The first specifies the
$\mathbf{k}$ (i.e., gives the ordinal corresponding to its position in
the list of k-points in `seedname.win`). The 2nd to 5th integers specify
$\mathbf{k+b}$. The 2nd integer, in particular, points to the k-point on
the list that is a periodic image of $\mathbf{k+b}$, and in particular
is the image that is actually mentioned in the list. The last three
integers specify the $\mathbf{G}$ vector, in reciprocal lattice units,
that brings the k-point specified by the second integer, and that thus
lives inside the first BZ zone, to the actual $\mathbf{k+b}$ that we
need.

Subsequent $\verb#num_bands#\times\verb#num_bands#$ lines of each block:
two real numbers per line. These are the real and imaginary parts,
respectively, of the actual scalar product $M_{mn}^{(\mathbf{k,b})}$ for
$m,n \in [1,\verb#num_bands#]$. The order of these elements is such that
the first index $m$ is fastest.

### `seedname.amn` file

INPUT.

The file `seedname.amn` contains the projection $A_{mn}^{(\mathbf{k})}$.

First line: a user comment, e.g., the date and time

Second line: 3 integers: `num_bands`, `num_kpts`, `num_wann`

Subsequently
$\verb#num_bands#\times\verb#num_wann#\times\verb#num_kpts#$ lines: 3
integers and 2 real numbers on each line. The first two integers are the
band index $m$ and the projection index $n$, respectively. The third
integer specifies the $\mathbf{k}$ by giving the ordinal corresponding
to its position in the list of $k$-points in `seedname.win`. The real
numbers are the real and imaginary parts, respectively, of the actual
$A_{mn}^{(\mathbf{k})}$.

### `seedname.dmn` file

INPUT.

The file `seedname.dmn` contains the data needed to construct
symmetry-adapted Wannier functions [@sakuma-prb13]. Required if
`site_symmetry = .true.`

First line: a user comment, e.g., the date and time

Second line: 4 integers: `num_bands`, `nsymmetry`, `nkptirr`,
`num_kpts`.\
`nsymmetry`: the number of symmetry operations\
`nkptirr`: the number of irreducible k-points

Blank line

`num_kpts` integers: Mapping between full k- and irreducible k-points.
Each k-point is related to some k-point in the irreducible BZ. The
information of this mapping is written. Each entry corresponds to a
k-point in the full BZ, in the order in which they appear in the k-point
list in `seedname.win` file. The (integer) value of each entry is the
k-point index in the IBZ to which the k-point maps. The number of unique
values is equal to the number of k-points in the IBZ. The data is
written 10 values per line.

Blank line

`nkptirr` integers: List of irreducible k-points. Each entry corresponds
to a k-point of the IBZ. The (integer) value of each entry is the
k-point index corresponding to the k-point list in `seedname.win` file.
The values should be between 1 and `num_kpts`. The data is written 10
values per line.

Blank line

`nkptirr` blocks of `nsymmetry` integer data (each block separated by a
blank line): List of k-points obtained by acting the symmetry operations
on the irreducible k-points. The data is written 10 values per line.

Blank line

$\verb#nsymmetry# \times \verb#nkptirr#$ blocks of data:\
The information of $D$ matrix in Eq. (15) of Ref. [@sakuma-prb13]. Each
block contains $\verb#num_wann# \times \verb#num_wann#$ lines and is
separated by a blank line. The data are stored in
`d_matrix_wann(m,n,isym,ikirr)` with
$\verb#m#, \verb#n# \in [1,\verb#num_wann#]$,
$\verb#isym# \in [1,\verb#nsymmetry#]$, and
$\verb#ikirr# \in [1,\verb#nkptirr#]$. The order of the elements is such
that left indices run faster than right indices (`m`: fastest, `ikirr`:
slowest).

Blank line

$\verb#nsymmetry# \times \verb#nkptirr#$ blocks of data:\
The information of $\tilde d$ matrix in Eq. (17) of
Ref. [@sakuma-prb13]. Each block contains
$\verb#num_bands# \times \verb#num_bands#$ lines and is separated by a
blank line. The data are stored in `d_matrix_band(m,n,isym,ikirr)` with
$\verb#m#, \verb#n# \in [1,\verb#num_bands#]$,
$\verb#isym# \in [1,\verb#nsymmetry#]$, and
$\verb#ikirr# \in [1,\verb#nkptirr#]$. The order of the elements is such
that left indices run faster than right indices (`m`: fastest, `ikirr`:
slowest).

### `seedname.eig` file

INPUT.

Required if any of `disentanglement`, `plot_bands`, `plot_fermi_surface`
or `write_hr` are `.true.`

The file `seedname.eig` contains the Kohn-Sham eigenvalues
$\varepsilon_{n\mathbf{k}}$ (in eV) at each point in the Monkhorst-Pack
mesh.

Each line consist of two integers and a real number. The first integer
is the band index, the second integer gives the ordinal corresponding to
the $k$-point in the list of $k$-points in `seedname.win`, and the real
number is the eigenvalue.

E.g.,

               1           1  -6.43858831271328
               2           1   19.3977795287297
               3           1   19.3977795287297
               4           1   19.3977795287298

### Interface with pwscf

Interfaces between `wannier90` and many ab-initio codes such as
[pwscf]{.smallcaps}, abinit (<http://www.abinit.org>), siesta
(<http://www.icmab.es/siesta/>), fleur, VASP and Wien2k
(<http://www.wien2k.at>) are available. Here we describe the seamless
interface between `wannier90` and [pwscf]{.smallcaps}, a plane-wave DFT
code that comes as part of the Quantum ESPRESSO package (see
<http://www.quantum-espresso.org>). You will need to download and
compile [pwscf]{.smallcaps} (i.e., the `pw.x` code) and the
post-processing interface `pw2wannier90.x`. Please refer to the
documentation that comes with the Quantum ESPRESSO distribution for
instructions.

1.  Run 'scf'/'nscf' calculation(s) with `pw`

2.  Run `wannier90` with `postproc_setup` = `.true.` to generate
    `seedname.nnkp`

3.  Run `pw2wannier90`. First it reads an input file, e.g.,
    `seedname.pw2wan`, which defines `prefix` and `outdir` for the
    underlying 'scf' calculation, as well as the name of the file
    `seedname.nnkp`, and does a consistency check between the direct and
    reciprocal lattice vectors read from `seedname.nnkp` and those
    defined in the files specified by `prefix`. `pw2wannier90` generates
    `seedname.mmn`, `seedname.amn` and `seedname.eig`. `seedname.dmn`
    and `seedname.sym` files are additionally created when
    `write_dmn = .true.` (see below).

4.  Run `wannier90` with `postproc_setup` = `.false.` to disentangle
    bands (if required), localise MLWF, and use MLWF for plotting,
    bandstructures, Fermi surfaces etc.

Examples of how the interface with [pwscf]{.smallcaps} works are given
in the `wannier90` Tutorial.

#### `seedname.pw2wan`

A number of keywords may be specified in the `pw2wannier90` input file:

-   `outdir` -- Location to write output files. Default is `` `./' ``

-   `prefix` -- Prefix for the [pwscf]{.smallcaps} calculation. Default
    is `` ` ' ``

-   `seedname` -- Seedname for the `wannier90` calculation. Default is
    `` `wannier' ``

-   `spin_component` -- Spin component. Takes values `` `up' ``,
    `` `down' `` or `` `none' `` (default).

-   `wan_mode` -- Either `` `standalone' `` (default) or `` `library' ``

-   `write_unk` -- Set to `.true.` to write the periodic part of the
    Bloch functions for plotting in `wannier90`. Default is `.false.`

-   `reduce_unk` -- Set to `.true.` to reduce file-size (and resolution)
    of Bloch functions by a factor of 8. Default is `.false.` (only
    relevant if `write_unk=.true.`)[^5]

-   `wvfn_formatted` -- Set to `.true.` to write formatted
    wavefunctions. Default is `.false.` (only relevant if
    `write_unk=.true.`)

-   `write_amn` -- Set to `.false.` if $A_{mn}^{(\mathbf{k})}$ not
    required. Default is `.true.`

-   `write_mmn` -- Set to `.false.` if $M_{mn}^{(\mathbf{k,b})}$ not
    required. Default is `.true.`

-   `write_spn` -- Set to `.true.` to write out the matrix elements of
    $S$ between Bloch states (non-collinear spin calculation only).
    Default is `.false.`

-   `spn_formatted` -- Set to `.true.` to write spn data as a formatted
    file. Default is `.false.` (only relevant if `write_spn=.true.`)

-   `write_uHu` -- Set to `.true.` to write out the matrix elements
    $$\langle u_{n{\bf k}+{\bf b}_1}\vert
    H_{\bf k}\vert u_{m{\bf k}+{\bf b}_2}\rangle.$$ Default is `.false.`

-   `uHu_formatted` -- Set to `.true.` to write uHu data as a formatted
    file. Default is `.false.` (only relevant if `write_uHu=.true.`)

-   `write_uIu` -- Set to `.true.` to write out the matrix elements of
    $$\langle  u_{n{\bf k}+{\bf b}_1}\vert
    u_{m{\bf k}+{\bf b}_2}\rangle.$$ Default is `.false.`

-   `uIu_formatted` -- Set to `.true.` to write uIu data as a formatted
    file. Default is `.false.` (only relevant if `write_uIu=.true.`)

-   `write_unkg` -- Set to `.true.` to write the first few Fourier
    components of the periodic parts of the Bloch functions.

-   `write_dmn` -- Set to `.true.` to construct symmetry-adapted Wannier
    functions. Default is `.false.`

-   `read_sym` -- Set to `.true.` to customize symmetry operations to be
    used in symmetry-adapted mode. When `read_sym = .true.`, an
    additional input `seedname.sym` is required. Default is `.false.`
    (only relevant if `write_dmn=.true.`).

-   `atom_proj` -- Set to `.true.` to use pseudo-atomic orbitals for
    computing `amn`. Default is `.false.`.

-   `atom_proj_exclude` -- A list of integers specifying the indices of
    pseudo-atomic projectors to be excluded from computing `amn`. Used
    only when `atom_proj = .true.`. No default.

-   `atom_proj_ext` -- Set to `.true.` to use external pseudo-atomic
    orbitals for computing `amn`, will read data files from directory
    `atom_proj_dir`. Used only when `atom_proj = .true.`. Default is
    `.false.`.

-   `atom_proj_dir` -- A string specifying the directory for external
    pseudo-atomic projectors. Used only when `atom_proj = .true.` and
    `atom_proj_ext = .true.`. No default.

For examples of use, refer to the `wannier90` Tutorial.

#### `seedname.sym`

If `read_sym = .true.`, then this additional input file is required for
`pw2wannier90.x`\
if `read_sym = .false.`, then this file is written by `pw2wannier90.x`
(only for reference -- it is not used in subsequent calculations)

The file `seedname.sym` contains the information of symmetry operations
used to create symmetry-adapted Wannier functions. If
`read_sym = .false.` (default), `pw2wannier90.x` uses the full symmetry
recognized by `pw.x`. If `read_sym = .true.`, you can specify symmetry
operations to be used in symmetry-adapted mode.

First line: an integer: `nsymmetry` (number of symmetry operations)

Second line: blank

Then: `nsymmetry` blocks of data. Each block (separated by a blank line)
consists of four lines. The order of the data in each block is as
follows:

      R(1,1)   R(2,1)   R(3,1)
      R(1,2)   R(2,2)   R(3,2)
      R(1,3)   R(2,3)   R(3,3)
       t(1)     t(2)     t(3)   

Here, $R$ is the rotational part of symmetry operations ($3\times3$
matrix), and $\bf t$ is the fractional translation in the unit of
"`alat`" (refer the definition of "`alat`" to the manual of
[pwscf]{.smallcaps}). Both data are given in Cartesian coordinates. The
symmetry operations act on a point $\bf r$ as ${\bf r} R - {\bf t}$.

## `wannier90` as a library {#ch:wann-lib}

This is a description of the interface between any external program and
the wannier code. There are two subroutines: `wannier_setup` and
`wannier_run`. Calling `wannier_setup` will return information required
to construct the $M_{mn}^{(\mathbf{k,b})}$ overlaps
(Ref. [@marzari-prb97], Eq. (25)) and
$A_{mn}^{(\mathbf{k})}=\left\langle
  \psi_{m\mathbf{k}}|g_{n}\right\rangle$ projections
(Ref. [@marzari-prb97], Eq. (62); Ref. [@souza-prb01], Eq. (22)). Once
the overlaps and projection have been computed, calling `wannier_run`
activates the minimisation and plotting routines in `wannier90`.

**IMPORTANT NOTE:** the library mode ONLY works in serial. Please call
it from a serial code, or if compiled in parallel, make sure to run it
from a single MPI process.

You can find a minimal example of how the library mode can be used among
the tests, in the file `test-suite/library-mode-test/test_library.F90`
in the Wannier90 git repository.

### Subroutines

#### `wannier_setup`

**`wannier_setup(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice,`\
`              kpt_latt,num_bands_tot,num_atoms,atom_symbols,atoms_cart,`\
`              gamma_only,spinors,nntot,nnlist,nncell,num_bands,num_wann,proj_site,`\
`              proj_l,proj_m,proj_radial,proj_z,proj_x,proj_zona,`\
`              exclude_bands,proj_s,proj_s_qaxis)`**

Conditions:

-   $\verb#num_kpts# = \verb#mp_grid(1)# \times \verb#mp_grid(2)#
    \times \verb#mp_grid(3)#$.

-   $\verb#num_nnmax# = 12$

This subroutine returns the information required to determine the
required overlap elements $M_{mn}^{(\mathbf{k,b})}$ and projections
$A_{mn}^{(\mathbf{k})}$, i.e., `M_matrix` and `A_matrix`, described in
Section [6.1.2](#wannier_run){reference-type="ref"
reference="wannier_run"}.

For the avoidance of doubt, `real_lattice(1,2)` is the $y-$component of
the first lattice vector $\mathbf{A}_{1}$, etc.

The list of nearest neighbours of a particular k-point `nkp` is given by
`nnlist(nkp,1:nntot)`.

Additionally, the parameter `shell_list` may be specified in the
`wannier90` input file.

#### `wannier_run`

**`wannier_run(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice,`\
`            kpt_latt,num_bands,num_wann,nntot,num_atoms,atom_symbols,`\
`            atoms_cart,gamma_only,M_matrix_orig,A_matrix,eigenvalues,`\
`            U_matrix,U_matrix_opt,lwindow,wann_centres,wann_spreads,`\
`            spread`)**

-   `character(len=*), intent(in) :: seed_name`\
    The seedname of the current calculation.

-   `integer, dimension(3), intent(in) :: mp_grid`\
    The dimensions of the Monkhorst-Pack k-point grid.

-   `integer, intent(in) :: num_kpts`\
    The number of k-points on the Monkhorst-Pack grid.

-   `real(kind=dp), dimension(3,3),` ` intent(in) :: real_lattice`\
    The lattice vectors in Cartesian co-ordinates in units of Angstrom.

-   `real(kind=dp), dimension(3,3), intent(in) :: recip_lattice`\
    The reciprical lattice vectors in Cartesian co-ordinates in units of
    inverse Angstrom.

-   `real(kind=dp), dimension(3,num_kpts),` ` intent(in) :: kpt_latt`\
    The positions of the k-points in fractional co-ordinates relative to
    the reciprocal lattice vectors.

-   `integer, intent(in) :: num_bands`\
    The total number of bands to be processed.

-   `integer, intent(in) :: num_wann`\
    The number of MLWF to be extracted.

-   `integer, intent(in) :: nntot`\
    The number of nearest neighbours for each k-point.

-   `integer, intent(in) :: num_atoms`\
    The total number of atoms in the system.

-   `character(len=20), dimension(num_atoms),`
    ` intent(in) :: atom_symbols`\
    The elemental symbols of the atoms.

-   `real(kind=dp), dimension(3,num_atoms),` `intent(in) :: atoms_cart`\
    The positions of the atoms in Cartesian co-ordinates in Angstrom.

-   `logical, intent(in) :: gamma_only`\
    Set to `.true.` if the underlying electronic structure calculation
    has been performed with only $\Gamma$-point sampling and, hence, if
    the Bloch eigenstates that are used to construct
    $A_{mn}^{(\mathbf{k})}$ and $M_{mn}^{\mathbf{(k,b)}}$ are real.

-   `complex(kind=dp),`
    ` dimension(num_bands,num_bands,nntot,num_kpts),`\
    `                  intent(in) :: M_matrix`\
    The matrices of overlaps between neighbouring periodic parts of the
    Bloch eigenstates at each k-point, $M_{mn}^{(\mathbf{(k,b)})}$
    (Ref. [@marzari-prb97], Eq. (25)).

-   `complex(kind=dp), dimension(num_bands,num_wann,num_kpts),`\
    `                  intent(in) :: A_matrix`\
    The matrices describing the projection of `num_wann` trial orbitals
    on `num_bands` Bloch states at each k-point, $A_{mn}^{(\mathbf{k})}$
    (Ref. [@marzari-prb97], Eq. (62); Ref. [@souza-prb01], Eq. (22)).

-   `real(kind=dp), dimension(num_bands,num_kpts),`
    `intent(in) :: eigenvalues`\
    The eigenvalues $\varepsilon_{n\mathbf{k}}$ corresponding to the
    eigenstates, in eV.

-   `complex(kind=dp), dimension(num_wann,num_wann,num_kpts),`\
    `                  intent(out) :: U_matrix`\
    The unitary matrices at each k-point (Ref. [@marzari-prb97],
    Eq. (59))

-   `complex(kind=dp), dimension(num_bands,num_wann,num_kpts),`\
    `               optional, intent(out) :: U_matrix_opt`\
    The unitary matrices that describe the optimal sub-space at each
    k-point (see Ref. [@souza-prb01], Section IIIa). The array is packed
    (see below)

-   `logical, dimension(num_bands,num_kpts), optional, intent(out) :: lwindow`\
    The element `lwindow(nband,nkpt)` is `.true.` if the band `nband`
    lies within the outer energy window at kpoint `nkpt`.

-   `real(kind=dp), dimension(3,num_wann), optional, intent(out) :: wann_centres`\
    The centres of the MLWF in Cartesian co-ordinates in Angstrom.

-   `real(kind=dp), dimension(num_wann), optional, intent(out) :: wann_spreads`\
    The spread of each MLWF in Å$^{2}$.

-   `real(kind=dp), dimension(3), optional, intent(out) ::` `spread`\
    The values of $\Omega$, $\Omega_{\mathrm{I}}$ and $\tilde{\Omega}$
    (Ref. [@marzari-prb97], Eq. (13)).

Conditions:

-   $\verb#num_wann# \le \verb#num_bands#$

-   $\verb#num_kpts# = \verb#mp_grid(1)# \times \verb#mp_grid(2)#
    \times \verb#mp_grid(3)#$.

If $\verb#num_bands# = \verb#num_wann#$ then `U_matrix_opt` is the
identity matrix and `lwindow=.true.`

For the avoidance of doubt, `real_lattice(1,2)` is the $y-$component of
the first lattice vector $\mathbf{A}_{1}$, etc.

$$\begin{aligned}
\verb#M_matrix(m,n,nn,nkp)# & = & \left\langle u_{m\mathbf{k}} |
u_{n\mathbf{k+b}}\right\rangle\\
\verb#A_matrix(m,n,nkp)# & = &
\left\langle \psi_{m\mathbf{k}}|g_{n}\right\rangle\\
\verb#eigenvalues(n,nkp)# &=& \varepsilon_{n\mathbf{k}}
\end{aligned}$$ where $$\begin{aligned}
\mathbf{k} &=&\verb#kpt_latt(1:3,nkp)#\\
\mathbf{k+b}&=& \verb#kpt_latt(1:3,nnlist(nkp,nn))# +
\verb#nncell(1:3,nkp,nn)# 
\end{aligned}$$ and $\left\{|g_{n}\rangle\right\}$ are a set of initial
trial orbitals. These are typically atom or bond-centred Gaussians that
are modulated by appropriate spherical harmonics.

Additional parameters should be specified in the `wannier90` input file.

## Transport Calculations with `wannier90`  {#ch:transport}

By setting $\verb#transport#=\verb#TRUE#$, `wannier90` will calculate
the quantum conductance and density of states of a one-dimensional
system. The results will be written to files `seedname_qc.dat` and
`seedname_dos.dat`, respectively.

The system for which transport properties are calculated is determined
by the keyword `transport_mode`.

### `transport_mode = bulk`

Quantum conductance and density of states are calculated for a perfectly
periodic one-dimensional conductor. If
$\verb#tran_read_ht#=\verb#FALSE#$ the transport properties are
calculated using the Hamiltonian in the Wannier function basis of the
system found by `wannier90`. Setting $\verb#tran_read_ht#=\verb#TRUE#$
allows the user to provide an external Hamiltonian matrix file
`seedname_htB.dat`, from which the properties are found. See
Section [2.9](#sec:post-p){reference-type="ref" reference="sec:post-p"}
for more details of the keywords required for such calculations.

### `transport_mode = lcr`

Quantum conductance and density of states are calculated for a system
where semi-infinite, left and right leads are connected through a
central conductor region. This is known as the *lcr* system. Details of
the method is described in Ref. [@nardelli-prb99].

In `wannier90` two options exist for performing such calculations:

-   If $\verb#tran_read_ht#=\verb#TRUE#$ the external Hamiltonian files
    `seedname_htL.dat, seedname_htLC.dat, seedname_htC.dat, seedname_htCR.dat, seedname_htR.dat`
    are read and used to compute the transport properties.

-   If $\verb#tran_read_ht#=\verb#FALSE#$, then the transport
    calculation is performed automatically using the Wannier functions
    as a basis and the 2c2 geometry described in
    Section [7.3](#sec:2c2){reference-type="ref" reference="sec:2c2"}.

### Automated lcr Transport Calculations: The 2c2 Geometry {#sec:2c2}

Calculations using the 2c2 geometry provide a method to calculate the
transport properties of an lcr system from a single
`wannier90` calculation. The Hamiltonian matrices which the five
external files provide in the $\verb#tran_read_ht#=\verb#TRUE#$ case are
instead built from the Wannier function basis directly. As such, strict
rules apply to the system geometry, which is shown in
Figure [7.1](#fig:2c2){reference-type="ref" reference="fig:2c2"}. These
rules are as follows:

-   Left and right leads must be identical and periodic.

-   Supercell must contain two principal layers (PLs) of lead on the
    left, a central conductor region and two principal layers of lead on
    the right.

-   The conductor region must contain enough lead such that the disorder
    does not affect the principal layers of lead either side.

-   A single **k**-point (Gamma) must be used.

![Schematic illustration of the supercell required for 2c2 lcr
calculations, showing where each of the Hamiltonian matrices are derived
from. Four principal layers (PLs) are required plus the conductor
region.](lcr_2c2.pdf){#fig:2c2 height="4cm"}

In order to build the Hamiltonians, Wannier functions are first sorted
according to position and then type if a number of Wannier functions
exist with a similar centre (eg. *d*-orbital type Wannier functions
centred on a Cu atom). Next, consistent parities of Wannier function are
enforced. To distingiush between different types of Wannier function and
assertain relative parities, a signature of each Wannier function is
computed. The signature is formed of 20 integrals which have different
spatial dependence. They are given by:

$$I=\frac{1}{V}\int_V g(\mathbf{r})w(\mathbf{r})d\mathbf{r}
\label{eq:sig_ints}$$

where $V$ is the volume of the cell, $w(\mathbf{r})$ is the Wannier
function and $g(\mathbf{r})$ are the set of functions:

$$\begin{aligned}
g(\mathbf{r})=&\left\lbrace1,\sin\left(\frac{2\pi (x-x_c)}{L_x}\right),
											 \sin\left(\frac{2\pi (y-y_c)}{L_y}\right),
											 \sin\left(\frac{2\pi (z-z_c)}{L_z}\right),
											 \sin\left(\frac{2\pi (x-x_c)}{L_x}\right)
											 \sin\left(\frac{2\pi (y-y_c)}{L_y}\right),\right.\nonumber \\
										   &\left.\sin\left(\frac{2\pi (x-x_c)}{L_x}\right)
											 \sin\left(\frac{2\pi (z-z_c)}{L_z}\right),
											 ... \right\rbrace
\label{eq:g(r)}
\end{aligned}$$ upto third order in powers of sines. Here, the supercell
has dimension $(L_x,L_y,L_z)$ and the Wannier function has centre
$\mathbf{r}_c=(x_c,y_c,z_c)$. Each of these integrals may be written as
linear combinations of the following sums:

$$S_n(\mathbf{G})=\displaystyle{e^{i\mathbf{G.r}_{c}}\sum_{m}U_{mn}\tilde{u}_{m\Gamma}^{*}(\mathbf{G})}$$

where $n$ and $m$ are the Wannier function and band indexes,
$\mathbf{G}$ is a G-vector, $U_{mn}$ is the unitary matrix that
transforms from the Bloch reopresentation of the system to the
maximally-localised Wannier function basis and
$\tilde{u}_{m\Gamma}^{*}(\mathbf{G})$ are the conjugates of the Fourier
transforms of the periodic parts of the Bloch states at the $\Gamma\!$
-point. The complete set of $\tilde{u}_{m\mathbf{k}}(\mathbf{G})$ are
often outputted by plane-wave DFT codes. However, to calculate the 20
signature integrals, only 32 specific
$\tilde{u}_{m\mathbf{k}}(\mathbf{G})$ are required. These are found in
an additional file (`seedname.unkg`) that should be provided by the
interface between the DFT code and `wannier90` . A detailed description
of this file may be found in
Section [8.32](#sec:files_unkg){reference-type="ref"
reference="sec:files_unkg"}.

Additionally, the following keywords are also required in the input
file:

-   `tran_num_ll` : The number of Wannier functions in a principal
    layer.

-   `tran_num_cell_ll` : The number of unit cells in one principal layer
    of lead

A further parameter related to these calculations is
`tran_group_threshold`.

Examples of how 2c2 calculations are preformed can be found in the
`wannier90` Tutorial.

## Files

### `seedname.win`

INPUT. The master input file; contains the specification of the system
and any parameters for the run. For a description of input parameters,
see Chapter [2](#chap:parameters){reference-type="ref"
reference="chap:parameters"}; for examples, see
Section [10.1](#winfile){reference-type="ref" reference="winfile"} and
the `wannier90` Tutorial.

#### Units

The following are the dimensional quantities that are specified in the
master input file:

-   Direct lattice vectors

-   Positions (of atomic or projection) centres in real space

-   Energy windows

-   Positions of k-points in reciprocal space

-   Convergence thresholds for the minimisation of $\Omega$

-   `zona` (see Section [3.1](#sec:proj){reference-type="ref"
    reference="sec:proj"})

-   `wannier_plot_cube`: cut-off radius for plotting WF in Gaussian cube
    format

Notes:

-   The units (either `ang` (default) or `bohr`) in which the lattice
    vectors, atomic positions or projection centres are given can be set
    in the first line of the blocks `unit_cell_cart`, `atoms_cart` and
    `projections`, respectively, in `seedname.win`.

-   Energy is always in eV.

-   Convergence thresholds are always in Å$^{2}$

-   Positions of k-points are always in crystallographic coordinates
    relative to the reciprocal lattice vectors.

-   `zona` is always in reciprocal Angstrom (Å$^{-1}$)

-   The keyword `length_unit` may be set to `ang` (default) or `bohr`,
    in order to set the units in which the quantities in the output file
    `seedname.wout` are written.

-   `wannier_plot_radius` is in Angstrom

The reciprocal lattice vectors
$\{\mathbf{B}_{1},\mathbf{B}_{2},\mathbf{B}_{3}\}$ are defined in terms
of the direct lattice vectors
$\{\mathbf{A}_{1},\mathbf{A}_{2},\mathbf{A}_{3}\}$ by the equation

$$\mathbf{B}_{1} = \frac{2\pi}{\Omega}\mathbf{A}_{2}\times\mathbf{A}_{3}
\ \ \ \mathrm{etc.},$$

where the cell volume is
$V=\mathbf{A}_{1}\cdot(\mathbf{A}_{2}\times\mathbf{A}_{3})$.

### `seedname.mmn`

INPUT. Written by the underlying electronic structure code. See
Chapter [5](#ch:wann-pp){reference-type="ref" reference="ch:wann-pp"}
for details.

### `seedname.amn`

INPUT. Written by the underlying electronic structure code. See
Chapter [5](#ch:wann-pp){reference-type="ref" reference="ch:wann-pp"}
for details.

### `seedname.dmn`

INPUT. Read if `site_symmetry = .true.` (symmetry-adapted mode). Written
by the underlying electronic structure code. See
Chapter [5](#ch:wann-pp){reference-type="ref" reference="ch:wann-pp"}
for details.

### `seedname.eig`

INPUT. Written by the underlying electronic structure code. See
Chapter [5](#ch:wann-pp){reference-type="ref" reference="ch:wann-pp"}
for details.

### `seedname.nnkp` {#sec:old-nnkp}

OUTPUT. Written by `wannier90` when `postproc_setup=.TRUE.` (or,
alternatively, when `wannier90` is run with the `-pp` command-line
option). See Chapter [5](#ch:wann-pp){reference-type="ref"
reference="ch:wann-pp"} for details.

### `seedname.wout`

OUTPUT. The master output file. Here we give a description of the main
features of the output. The verbosity of the output is controlled by the
input parameter `iprint`. The higher the value, the more detail is given
in the output file. The default value is 1, which prints minimal
information.

#### Header

The header provides some basic information about `wannier90`, the
authors, the code version and release, and the execution time of the
current run. The header looks like the following different (the string
might slightly change across different versions):


                 +---------------------------------------------------+
                 |                                                   |
                 |                   WANNIER90                       |
                 |                                                   |
                 +---------------------------------------------------+
                 |                                                   |
                 |        Welcome to the Maximally-Localized         |
                 |        Generalized Wannier Functions code         |
                 |            http://www.wannier.org                 |
                 |                                                   |
                 |  Wannier90 Developer Group:                       |
                 |    Giovanni Pizzi    (EPFL)                       |
                 |    Valerio Vitale    (Cambridge)                  |
                 |    David Vanderbilt  (Rutgers University)         |
                 |    Nicola Marzari    (EPFL)                       |
                 |    Ivo Souza         (Universidad del Pais Vasco) |
                 |    Arash A. Mostofi  (Imperial College London)    |
                 |    Jonathan R. Yates (University of Oxford)       |
                 |                                                   |
                 |  For the full list of Wannier90 3.x authors,      |
                 |  please check the code documentation and the      |
                 |  README on the GitHub page of the code            |
                 |                                                   |
                 |                                                   |
                 |  Please cite                                      |
                                           .
                                           .
                 |                                                   |
                 +---------------------------------------------------+
                 |    Execution started on 18Dec2018 at 18:39:42     |
                 +---------------------------------------------------+

#### System information

This part of the output file presents information that `wannier90` has
read or inferred from the master input file `seedname.win`. This
includes real and reciprocal lattice vectors, atomic positions,
k-points, parameters for job control, disentanglement, localisation and
plotting.

                                        ------
                                        SYSTEM
                                        ------
     
                                  Lattice Vectors (Ang)
                        a_1     3.938486   0.000000   0.000000
                        a_2     0.000000   3.938486   0.000000
                        a_3     0.000000   0.000000   3.938486
     
                       Unit Cell Volume:      61.09251  (Ang^3)
     
                            Reciprocal-Space Vectors (Ang^-1)
                        b_1     1.595330   0.000000   0.000000
                        b_2     0.000000   1.595330   0.000000
                        b_3     0.000000   0.000000   1.595330
      
     *----------------------------------------------------------------------------*
     |   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |
     +----------------------------------------------------------------------------+
     | Ba   1   0.00000   0.00000   0.00000   |    0.00000   0.00000   0.00000    |
     | Ti   1   0.50000   0.50000   0.50000   |    1.96924   1.96924   1.96924    |
                                              .
                                              . 
     *----------------------------------------------------------------------------*
      
                                    ------------
                                    K-POINT GRID
                                    ------------
      
                 Grid size =  4 x  4 x  4      Total points =   64
      
     *---------------------------------- MAIN ------------------------------------*
     |  Number of Wannier Functions               :                 9             |
     |  Number of input Bloch states              :                 9             |
     |  Output verbosity (1=low, 5=high)          :                 1             |
     |  Length Unit                               :               Ang             |
     |  Post-processing setup (write *.nnkp)      :                 F             |
                                                  .
                                                  .
     *----------------------------------------------------------------------------*

#### Nearest-neighbour k-points

This part of the output files provides information on the
$\mathrm{b}$-vectors and weights chosen to satisfy the condition of
Eq. [\[eq:B1\]](#eq:B1){reference-type="ref" reference="eq:B1"}.

     *---------------------------------- K-MESH ----------------------------------*
     +----------------------------------------------------------------------------+
     |                    Distance to Nearest-Neighbour Shells                    |
     |                    ------------------------------------                    |
     |          Shell             Distance (Ang^-1)          Multiplicity         |
     |          -----             -----------------          ------------         |
     |             1                   0.398833                      6            |
     |             2                   0.564034                     12            |
                                           .
                                           .
     +----------------------------------------------------------------------------+
     | The b-vectors are chosen automatically                                     |
     | The following shells are used:   1                                         |
     +----------------------------------------------------------------------------+
     |                        Shell   # Nearest-Neighbours                        |
     |                        -----   --------------------                        |
     |                          1               6                                 |
     +----------------------------------------------------------------------------+
     | Completeness relation is fully satisfied [Eq. (B1), PRB 56, 12847 (1997)]  |
     +----------------------------------------------------------------------------+

#### Disentanglement

Then (if required) comes the part where $\Omega_{\mathrm{I}}$ is
minimised to disentangle the optimally-connected subspace of states for
the localisation procedure in the next step.

First, a summary of the energy windows that are being used is given:

     *------------------------------- DISENTANGLE --------------------------------*
     +----------------------------------------------------------------------------+
     |                              Energy  Windows                               |
     |                              ---------------                               |
     |                   Outer:    2.81739  to   38.00000  (eV)                   |
     |                   Inner:    2.81739  to   13.00000  (eV)                   |
     +----------------------------------------------------------------------------+

Then, each step of the iterative minimisation of $\Omega_{\mathrm{I}}$
is reported.

                       Extraction of optimally-connected subspace                  
                       ------------------------------------------                  
     +---------------------------------------------------------------------+<-- DIS
     |  Iter     Omega_I(i-1)      Omega_I(i)      Delta (frac.)    Time   |<-- DIS
     +---------------------------------------------------------------------+<-- DIS
           1       3.82493590       3.66268867       4.430E-02      0.36    <-- DIS
           2       3.66268867       3.66268867       6.911E-15      0.37    <-- DIS
                                           .
                                           .
                                       
                 <<<      Delta < 1.000E-10  over  3 iterations     >>>
                 <<< Disentanglement convergence criteria satisfied >>>

            Final Omega_I     3.66268867 (Ang^2)

     +----------------------------------------------------------------------------+

The first column gives the iteration number. For a description of the
minimisation procedure and expressions for $\Omega_{\mathrm{I}}^{(i)}$,
see the original paper [@souza-prb01]. The procedure is considered to be
converged when the fractional difference between
$\Omega_{\mathrm{I}}^{(i)}$ and $\Omega_{\mathrm{I}}^{(i-1)}$ is less
than `dis_conv_tol` over ` dis_conv_window` iterations. The final column
gives a running account of the wall time (in seconds) so far. Note that
at the end of each line of output, there are the characters "`<– DIS`".
This enables fast searching of the output using, for example, the Unix
command `grep`:

`my_shell> grep DIS wannier.wout | less`

#### Wannierisation {#sec:files-wannierisation}

The next part of the output file provides information on the
minimisation of $\widetilde{\Omega}$. At each iteration, the centre and
spread of each WF is reported.

    *------------------------------- WANNIERISE ---------------------------------*
     +--------------------------------------------------------------------+<-- CONV
     | Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV
     +--------------------------------------------------------------------+<-- CONV
     
     ------------------------------------------------------------------------------
     Initial State
      WF centre and spread    1  (  0.000000,  1.969243,  1.969243 )     1.52435832
      WF centre and spread    2  (  0.000000,  1.969243,  1.969243 )     1.16120620
                                          .
                                          .
          0     0.126E+02     0.0000000000       12.6297685260       0.29  <-- CONV
            O_D=      0.0000000 O_OD=      0.1491718 O_TOT=     12.6297685 <-- SPRD
     ------------------------------------------------------------------------------
     Cycle:      1
      WF centre and spread    1  (  0.000000,  1.969243,  1.969243 )     1.52414024
      WF centre and spread    2  (  0.000000,  1.969243,  1.969243 )     1.16059775
                                          .
                                          .
      Sum of centres and spreads ( 11.815458, 11.815458, 11.815458 )    12.62663472
     
          1    -0.313E-02     0.0697660962       12.6266347170       0.34  <-- CONV
            O_D=      0.0000000 O_OD=      0.1460380 O_TOT=     12.6266347 <-- SPRD
     Delta: O_D= -0.4530841E-18 O_OD= -0.3133809E-02 O_TOT= -0.3133809E-02 <-- DLTA
     ------------------------------------------------------------------------------
     Cycle:      2
      WF centre and spread    1  (  0.000000,  1.969243,  1.969243 )     1.52414866
      WF centre and spread    2  (  0.000000,  1.969243,  1.969243 )     1.16052405
                                          .
                                          .
       Sum of centres and spreads ( 11.815458, 11.815458, 11.815458 )    12.62646411
     
          2    -0.171E-03     0.0188848262       12.6264641055       0.38  <-- CONV
            O_D=      0.0000000 O_OD=      0.1458674 O_TOT=     12.6264641 <-- SPRD
     Delta: O_D= -0.2847260E-18 O_OD= -0.1706115E-03 O_TOT= -0.1706115E-03 <-- DLTA
     ------------------------------------------------------------------------------
                                          .
                                          .
     ------------------------------------------------------------------------------
     Final State
      WF centre and spread    1  (  0.000000,  1.969243,  1.969243 )     1.52416618
      WF centre and spread    2  (  0.000000,  1.969243,  1.969243 )     1.16048545
                                          .
                                          .
      Sum of centres and spreads ( 11.815458, 11.815458, 11.815458 )    12.62645344
     
             Spreads (Ang^2)       Omega I      =    12.480596753
            ================       Omega D      =     0.000000000
                                   Omega OD     =     0.145856689
        Final Spread (Ang^2)       Omega Total  =    12.626453441
     ------------------------------------------------------------------------------

It looks quite complicated, but things look more simple if one uses
`grep`:

`my_shell> grep CONV wannier.wout`

gives

     +--------------------------------------------------------------------+<-- CONV
     | Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV
     +--------------------------------------------------------------------+<-- CONV
          0     0.126E+02     0.0000000000       12.6297685260       0.29  <-- CONV
          1    -0.313E-02     0.0697660962       12.6266347170       0.34  <-- CONV
                                                       .
                                                       .
         50     0.000E+00     0.0000000694       12.6264534413       2.14  <-- CONV

The first column is the iteration number, the second is the change in
$\Omega$ from the previous iteration, the third is the root-mean-squared
gradient of $\Omega$ with respect to variations in the unitary matrices
$\mathbf{U}^{(\mathbf{k})}$, and the last is the time taken (in
seconds). Depending on the input parameters used, the procedure either
runs for `num_iter` iterations, or a convergence criterion is applied on
$\Omega$. See Section [2.8](#sec:wann_params){reference-type="ref"
reference="sec:wann_params"} for details.

Similarly, the command

`my_shell> grep SPRD wannier.wout`

gives

            O_D=      0.0000000 O_OD=      0.1491718 O_TOT=     12.6297685 <-- SPRD
            O_D=      0.0000000 O_OD=      0.1460380 O_TOT=     12.6266347 <-- SPRD
                                                .
                                                .
            O_D=      0.0000000 O_OD=      0.1458567 O_TOT=     12.6264534 <-- SPRD         

which, for each iteration, reports the value of the diagonal and
off-diagonal parts of the non-gauge-invariant spread, as well as the
total spread, respectively. Recall from
Section [1](#sec:method){reference-type="ref" reference="sec:method"}
that
$\Omega = \Omega_{\mathrm{I}}+ \Omega_{\mathrm{D}} + \Omega_{\mathrm{OD}}$.

##### Wannierisation with selective localization and constrained centres

For full details of the selectively localised Wannier function (SLWF)
method, the reader is referred to Ref. [@Marianetti]. When using the
SLWF method, only a few things change in the output file and in general
the same principles described above will apply. In particular, when
minimising the spread with respect to the degrees of freedom of only a
subset of functions, it is not possible to cast the total spread
functional $\Omega$ as a sum of a gauge-invariant part and a
gauge-dependent part. Instead, one has
$\Omega^{'} = \Omega_{\mathrm{IOD}} + \Omega_{\mathrm{D}}$, where
$$\Omega^{'} = \sum_{n=1}^{J'<J} \left[\langle r^2 \rangle_n - \overline{\mathbf{r}}_{n}^{2}\right]$$
and
$$\Omega_{\mathrm{IOD}} = \sum_{n=1}^{J'<J} \left[\langle r^2_n \rangle- \sum_{\mathbf{R}} \vert\langle\mathbf{R}n\vert \mathbf{r} \vert n\mathbf{R}\rangle\vert^2 \right].$$
The total number of Wannier functions is $J$, whereas $J'$ is the number
functions to be selectively localized (so-called *objective WFs*). The
information on the number of functions which are going to be selectively
localized (`Number of Objective Wannier Functions`) is given in the
`MAIN` section of the output file:

     *---------------------------------- MAIN ------------------------------------*
     |  Number of Wannier Functions               :                 4             |
     |  Number of Objective Wannier Functions     :                 1             |
     |  Number of input Bloch states              :                 4             |

Whether or not the selective localization procedure has been switched on
is reported in the `WANNIERISE` section as

     |  Perform selective localization            :                 T             |

The next part of the output file provides information on the
minimisation of the modified spread functional:

     *------------------------------- WANNIERISE ---------------------------------*
     +--------------------------------------------------------------------+<-- CONV
     | Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV
     +--------------------------------------------------------------------+<-- CONV

     ------------------------------------------------------------------------------
     Initial State
      WF centre and spread    1  ( -0.857524,  0.857524,  0.857524 )     1.80463310
      WF centre and spread    2  (  0.857524, -0.857524,  0.857524 )     1.80463311
      WF centre and spread    3  (  0.857524,  0.857524, -0.857524 )     1.80463311
      WF centre and spread    4  ( -0.857524, -0.857524, -0.857524 )     1.80463311
      Sum of centres and spreads ( -0.000000, -0.000000,  0.000000 )     7.21853243

          0    -0.317E+01     0.0000000000       -3.1653368719       0.00  <-- CONV
           O_D=      0.0000000 O_IOD=     -3.1653369 O_TOT=     -3.1653369 <-- SPRD
     ------------------------------------------------------------------------------
     Cycle:      1
      WF centre and spread    1  ( -0.853260,  0.853260,  0.853260 )     1.70201498
      WF centre and spread    2  (  0.857352, -0.857352,  0.862454 )     1.84658331
      WF centre and spread    3  (  0.857352,  0.862454, -0.857352 )     1.84658331
      WF centre and spread    4  ( -0.862454, -0.857352, -0.857352 )     1.84658331
      Sum of centres and spreads ( -0.001010,  0.001010,  0.001010 )     7.24176492

          1    -0.884E-01     0.2093698260       -3.2536918930       0.00  <-- CONV
           O_IOD=     -3.2536919 O_D=      0.0000000 O_TOT=     -3.2536919 <-- SPRD
    Delta: O_IOD= -0.1245020E+00 O_D=  0.0000000E+00 O_TOT= -0.8835502E-01 <-- DLTA
     ------------------------------------------------------------------------------
                                          .
                                          .
     ------------------------------------------------------------------------------
     Final State
      WF centre and spread    1  ( -0.890189,  0.890189,  0.890189 )     1.42375495
      WF centre and spread    2  (  0.895973, -0.895973,  0.917426 )     2.14313664
      WF centre and spread    3  (  0.895973,  0.917426, -0.895973 )     2.14313664
      WF centre and spread    4  ( -0.917426, -0.895973, -0.895973 )     2.14313664
      Sum of centres and spreads ( -0.015669,  0.015669,  0.015669 )     7.85316486
     
             Spreads (Ang^2)       Omega IOD    =     1.423371553
            ================       Omega D      =     0.000383395
                                   Omega Rest   =     9.276919811
        Final Spread (Ang^2)       Omega Total  =     1.423754947
     ------------------------------------------------------------------------------

When comparing the output from an SLWF calculation with a standard
wannierisation (see
Sec. [8.7.5](#sec:files-wannierisation){reference-type="ref"
reference="sec:files-wannierisation"}), the only differences are in the
definition of the spread functional. Hence, during the minimization
`O_OD` is replaced by `O_IOD` and `O_TOT` now reflects the fact that the
new total spread functional is $\Omega^{'}$. The part on the final state
has one more item of information: the value of the difference between
the global spread functional and the new spread functional given by
`Omega Rest`
$$\Omega_{R} = \sum_{n=1}^{J-J'} \left[\langle r^2 \rangle_n - \overline{\mathbf{r}}_{n}^{2} \right]$$

If adding centre-constraints to the SLWFs, you will find the information
about the centres of the original projections and the desired centres in
the `SYSTEM` section

     *----------------------------------------------------------------------------*
     | Wannier#        Original Centres              Constrained centres          |
     +----------------------------------------------------------------------------+
     |    1     0.25000   0.25000   0.25000   |    0.00000   0.00000   0.00000    |
     *----------------------------------------------------------------------------*

As before one can check that the selective localization with constraints
is being used by looking at the `WANNIERISE` section:

     |  Perform selective localization            :                 T             |
     |  Use constrains in selective localization  :                 T             |
     |  Value of the Lagrange multiplier          :         0.100E+01             |
     *----------------------------------------------------------------------------*

which also gives the selected value for the Lagrange multiplier. The
output file for the minimisation section is modified as follows: both
`O_IOD` and `O_TOT` now take into account the factors coming from the
new term in the functional due to the constraints, which are implemented
by adding the following penalty functional to the spread functional,
$$\lambda_c \sum_{n=1}^{J'} \left(\overline{\mathbf{r}}_n - \mathbf{r}_{0n} \right)^2,$$
where $\mathbf{r}_{0n}$ is the desired centre for the $n^{\text{th}}$
Wannier function, see Ref. [@Marianetti] for details. The layout of the
output file at each iteration is unchanged.

          1    -0.884E-01     0.2093698260       -3.2536918930       0.00  <-- CONV

As regarding the final state, the only addition is the information on
the value of the penalty functional associated with the constraints
(`Penalty func`), which should be zero if the final centres of the
Wannier functions are at the target centres:

     Final State
      WF centre and spread    1  ( -1.412902,  1.412902,  1.412902 )     1.63408756
      WF centre and spread    2  (  1.239678, -1.239678,  1.074012 )     2.74801593
      WF centre and spread    3  (  1.239678,  1.074012, -1.239678 )     2.74801592
      WF centre and spread    4  ( -1.074012, -1.239678, -1.239678 )     2.74801592
      Sum of centres and spreads ( -0.007559,  0.007559,  0.007559 )     9.87813534

             Spreads (Ang^2)       Omega IOD_C   =    -4.261222001
            ================       Omega D       =     0.000000000
                                   Omega Rest    =     5.616913337
                                   Penalty func  =     0.000000000
        Final Spread (Ang^2)       Omega Total_C =    -4.261222001
     ------------------------------------------------------------------------------

#### Plotting

After WF have been localised, `wannier90` enters its plotting routines
(if required). For example, if you have specified an interpolated
bandstucture:

     *---------------------------------------------------------------------------*
     |                               PLOTTING                                    |
     *---------------------------------------------------------------------------*
      
     Calculating interpolated band-structure

#### Summary timings

At the very end of the run, a summary of the time taken for various
parts of the calculation is given. The level of detail is controlled by
the `timing_level` input parameter (set to 1 by default).

     *===========================================================================*
     |                             TIMING INFORMATION                            |
     *===========================================================================*
     |    Tag                                                Ncalls      Time (s)|
     |---------------------------------------------------------------------------|
     |kmesh: get                                        :         1         0.212|
     |overlap: read                                     :         1         0.060|
     |wann: main                                        :         1         1.860|
     |plot: main                                        :         1         0.168|
     *---------------------------------------------------------------------------*
     
     All done: wannier90 exiting

### `seedname.chk`

INPUT/OUTPUT. Information required to restart the calculation or enter
the plotting phase. If we have used disentanglement this file also
contains the rectangular matrices $\bf{U}^{{\rm dis}({\bf k})}$.

### `seedname.r2mn`

OUTPUT. Written if $\verb#write_r2mn#=\verb#true#$. The matrix elements
$\langle m|r^2|n\rangle$ (where $m$ and $n$ refer to MLWF)

### `seedname_band.dat`

OUTPUT. Written if `bands_plot=.TRUE.`; The raw data for the
interpolated band structure.

### `seedname_band.gnu`

OUTPUT. Written if `bands_plot=.TRUE.` and ` bands_plot_format=gnuplot`;
A `gnuplot` script to plot the interpolated band structure.

### `seedname_band.agr`

OUTPUT. Written if `bands_plot=.TRUE.` and ` bands_plot_format=xmgrace`;
A `grace` file to plot the interpolated band structure.

### `seedname_band.kpt`

OUTPUT. Written if `bands_plot=.TRUE.`; The k-points used for the
interpolated band structure, in units of the reciprocal lattice vectors.
This file can be used to generate a comparison band structure from a
first-principles code.

### `seedname.bxsf`

OUTPUT. Written if `fermi_surface_plot=.TRUE.`; A Fermi surface plot
file suitable for plotting with XCrySDen.

### `seedname_w.xsf`

OUTPUT. Written if `wannier_plot=.TRUE.` and
` wannier_plot_format=xcrysden`. Contains the ` w`$^{\mathrm{th}}$ WF in
real space in a format suitable for plotting with XCrySDen or VMD, for
example.

### `seedname_w.cube`

OUTPUT. Written if `wannier_plot=.TRUE.` and
` wannier_plot_format=cube`. Contains the ` w`$^{\mathrm{th}}$ WF in
real space in Gaussian cube format, suitable for plotting in XCrySDen,
VMD, gopenmol etc.

### `UNKp.s`

INPUT. Read if `wannier_plot`=`.TRUE.` and used to plot the MLWF. Read
if `transport_mode`=`lcr` and `tran_read_ht`=`.FALSE.` for use in
automated lcr transport calculations.

The periodic part of the Bloch states represented on a regular real
space grid, indexed by k-point `p` (from 1 to `num_kpts`) and spin `s`
('1' for 'up', '2' for 'down').

The name of the wavefunction file is assumed to have the form:

        write(wfnname,200) p,spin
    200 format ('UNK',i5.5,'.',i1)

The first line of each file should contain 5 integers: the number of
grid points in each direction (`ngx`, `ngy` and `ngz`), the k-point
number `ik` and the total number of bands `num_band` in the file. The
full file will be read by `wannier90` as:

    read(file_unit) ngx,ngy,ngz,ik,nbnd  
    do loop_b=1,num_bands
      read(file_unit) (r_wvfn(nx,loop_b),nx=1,ngx*ngy*ngz)
    end do

If `spinors`=`true` then `s`='NC', and the name of the wavefunction file
is assumed to have the form:

        write(wfnname,200) p
    200 format ('UNK',i5.5,'.NC')

and the file will be read by `wannier90` as:

    read(file_unit) ngx,ngy,ngz,ik,nbnd  
    do loop_b=1,num_bands
       read(file_unit) (r_wvfn_nc(nx,loop_b,1),nx=1,ngx*ngy*ngz) ! up-spinor
       read(file_unit) (r_wvfn_nc(nx,loop_b,2),nx=1,ngx*ngy*ngz) ! down-spinor
    end do  

All UNK files can be in formatted or unformatted style, this is
controlled by the logical keyword `wvfn_formatted`.

### `seedname_centres.xyz`

OUTPUT. Written if `write_xyz=.TRUE.`; xyz format atomic structure file
suitable for viewing with your favourite visualiser (`jmol`, `gopenmol`,
`vmd`, etc.).

### `seedname_hr.dat`

OUTPUT. Written if `write_hr=.TRUE.`. The first line gives the date and
time at which the file was created. The second line states the number of
Wannier functions `num_wann`. The third line gives the number of
Wigner-Seitz grid-points `nrpts`. The next block of `nrpts` integers
gives the degeneracy of each Wigner-Seitz grid point, with 15 entries
per line. Finally, the remaining `num_wann`$^2 \times$ `nrpts` lines
each contain, respectively, the components of the vector $\mathbf{R}$ in
terms of the lattice vectors $\{\mathbf{A}_{i}\}$, the indices $m$ and
$n$, and the real and imaginary parts of the Hamiltonian matrix element
$H_{mn}^{(\mathbf{R})}$ in the WF basis, e.g.,

     Created on 24May2007 at 23:32:09                            
            20
            17
        4   1   2    1    4    1    1    2    1    4    6    1    1   1   2
        1   2
        0   0  -2    1    1   -0.001013    0.000000
        0   0  -2    2    1    0.000270    0.000000
        0   0  -2    3    1   -0.000055    0.000000
        0   0  -2    4    1    0.000093    0.000000
        0   0  -2    5    1   -0.000055    0.000000
        .
        .
        .

### `seedname_r.dat`

OUTPUT. Written if $\verb#write_rmn#=\verb#true#$. The matrix elements
$\langle m\mathbf{0}|\mathbf{r}|n\mathbf{R}\rangle$ (where $n\mathbf{R}$
refers to MLWF $n$ in unit cell $\mathbf{R}$). The first line gives the
date and time at which the file was created. The second line states the
number of Wannier functions `num_wann`. The third line states the number
of $\mathbf{R}$ vectors `nrpts`. Similar to the case of the Hamiltonian
matrix above, the remaining `num_wann`$^2 \times$ `nrpts` lines each
contain, respectively, the components of the vector $\mathbf{R}$ in
terms of the lattice vectors $\{\mathbf{A}_{i}\}$, the indices $m$ and
$n$, and the real and imaginary parts of the position matrix element in
the WF basis.

### `seedname_tb.dat`

OUTPUT. Written if `write_tb=.TRUE.`. This file is essentially a
combination of `seedname_hr.dat` and `seedname_r.dat`, plus lattice
vectors. The first line gives the date and time at which the file was
created. The second to fourth lines are the lattice vectors in Angstrom
unit.

     written on 27Jan2020 at 18:08:42 
      -1.8050234585004898        0.0000000000000000        1.8050234585004898     
       0.0000000000000000        1.8050234585004898        1.8050234585004898     
      -1.8050234585004898        1.8050234585004898        0.0000000000000000 

The next part is the same as `seedname_hr.dat`. The fifth line states
the number of Wannier functions `num_wann`. The sixth line gives the
number of Wigner-Seitz grid-points `nrpts`. The next block of `nrpts`
integers gives the degeneracy of each Wigner-Seitz grid point, with 15
entries per line. Then, the next `num_wann`$^2 \times$ `nrpts` lines
each contain, respectively, the components of the vector $\mathbf{R}$ in
terms of the lattice vectors $\{\mathbf{A}_{i}\}$, the indices $m$ and
$n$, and the real and imaginary parts of the Hamiltonian matrix element
$H_{mn}^{(\mathbf{R})}$ in the WF basis, e.g.,

               7
              93
        4    6    2    2    2    1    2    2    1    1    2    6    2    2    2
        6    2    2    4    1    1    1    4    1    1    1    1    2    1    1
        1    2    2    1    1    2    4    2    1    2    1    1    1    1    2
        1    1    1    2    1    1    1    1    2    1    2    4    2    1    1
        2    2    1    1    1    2    1    1    1    1    4    1    1    1    4
        2    2    6    2    2    2    6    2    1    1    2    2    1    2    2
        2    6    4

       -3    1    1
        1    1    0.42351556E-02 -0.95722060E-07
        2    1    0.69481480E-07 -0.20318638E-06
        3    1    0.10966508E-06 -0.13983284E-06
        .
        .
        .

Finally, the last part is the same as `seedname_r.dat`. The
`num_wann`$^2 \times$ `nrpts` lines each contain, respectively, the
components of the vector $\mathbf{R}$ in terms of the lattice vectors
$\{\mathbf{A}_{i}\}$, the indices $m$ and $n$, and the real and
imaginary parts of the position matrix element in the WF basis (the
float numbers in columns 3 and 4 are the real and imaginary parts for
$\langle m\mathbf{0}|\mathbf{r}_x|n\mathbf{R}\rangle$, columns 5 and 6
for $\langle m\mathbf{0}|\mathbf{r}_y|n\mathbf{R}\rangle$, and columns 7
and 8 for $\langle m\mathbf{0}|\mathbf{r}_z|n\mathbf{R}\rangle$), e.g.

       -3    1    1
        1    1    0.32277552E-09  0.21174901E-08 -0.85436987E-09  0.26851510E-08  ...
        2    1   -0.18881883E-08  0.21786973E-08  0.31123076E-03  0.39228431E-08  ...
        3    1    0.31123242E-03 -0.35322230E-09  0.70867281E-09  0.10433480E-09  ...
        .
        .
        .

### `seedname.bvec`

OUTPUT. Written if $\verb#write_bvec#=\verb#true#$. This file contains
the matrix elements of bvector and their weights. The first line gives
the date and time at which the file was created. The second line states
the number of k-points and the total number of neighbours for each
k-point `nntot`. Then all the other lines contain the b-vector (x,y,z)
coordinate and weigths for each k-points and each of its neighbours.

### `seedname_wsvec.dat`

OUTPUT. Written if $\verb#write_hr#=\verb#true#$ or
$\verb#write_rmn#=\verb#true#$ or $\verb#write_tb#=\verb#true#$. The
first line gives the date and time at which the file was created and the
value of `use_ws_distance`. For each pair of Wannier functions
(identified by the components of the vector $\mathbf{R}$ separating
their unit cells and their indices) it gives: (i) the number of lattice
vectors of the periodic supercell $\mathbf{T}$ that bring the Wannier
function in $\mathbf{R}$ back in the Wigner-Seitz cell centred on the
other Wannier function and (ii) the set of superlattice vectors
$\mathbf{T}$ to make this transformation. These superlattice vectors
$\mathbf{T}$ should be added to the $\mathbf{R}$ vector to obtain the
correct centre of the Wannier function that underlies a given matrix
element (e.g. the Hamiltonian matrix elements in `seedname_hr.dat`) in
order to correctly interpolate in reciprocal space.

    ## written on 20Sep2016 at 18:12:37  with use_ws_distance=.true.
        0    0    0    1    1
        1
        0    0    0
        0    0    0    1    2
        1
        0    0    0
        0    0    0    1    3
        1
        0    0    0
        0    0    0    1    4
        1
        0    0    0
        0    0    0    1    5
        1
        0    0    0
        0    0    0    1    6
        2
        0   -1   -1
        1   -1   -1
        .
        .
        .  

### `seedname_qc.dat`

OUTPUT. Written if $\verb#transport#=\verb#.TRUE.#$. The first line
gives the date and time at which the file was created. In the subsequent
lines, the energy value in units of eV is written in the left column,
and the quantum conductance in units of $\frac{2e^2}{h}$
($\frac{e^2}{h}$ for a spin-polarized system) is written in the right
column.

     ## written on 14Dec2007 at 11:30:17
       -3.000000       8.999999
       -2.990000       8.999999
       -2.980000       8.999999
       -2.970000       8.999999
        .
        .
        .

### `seedname_dos.dat`

OUTPUT. Written if $\verb#transport#=\verb#.TRUE.#$. The first line
gives the date and time at which the file was created. In the subsequent
lines, the energy value in units of eV is written in the left column,
and the density of states in an arbitrary unit is written in the right
column.

     ## written on 14Dec2007 at 11:30:17
       -3.000000       6.801199
       -2.990000       6.717692
       -2.980000       6.640828
       -2.970000       6.569910
        .
        .
        .

### `seedname_htB.dat`

INPUT/OUTPUT. Read if $\verb#transport_mode#=\verb#bulk#$ and
$\verb#tran_read_ht#=\verb#.TRUE.#$. Written if
$\verb#tran_write_ht#=\verb#.TRUE.#$. The first line gives the date and
time at which the file was created. The second line gives `tran_num_bb`.
The subsequent lines contain `tran_num_bb`$\times$`tran_num_bb` $H_{mn}$
matrix, where the indices $m$ and $n$ span all `tran_num_bb` WFs located
at $0^{\mathrm{th}}$ principal layer. Then `tran_num_bb` is recorded
again in the new line followed by $H_{mn}$, where $m^{\mathrm{th}}$ WF
is at $0^{\mathrm{th}}$ principal layer and $n^{\mathrm{th}}$ at
$1^{\mathrm{st}}$ principal layer. The $H_{mn}$ matrix is written in
such a way that $m$ is the fastest varying index.

     written on 14Dec2007 at 11:30:17
       150
       -1.737841   -2.941054    0.052673   -0.032926    0.010738   -0.009515
        0.011737   -0.016325    0.051863   -0.170897   -2.170467    0.202254
        .
        .
        .
       -0.057064   -0.571967   -0.691431    0.015155   -0.007859    0.000474
       -0.000107   -0.001141   -0.002126    0.019188   -0.686423  -10.379876
       150
        0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
        0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
        .
        .
        .
        0.000000    0.000000    0.000000    0.000000    0.000000   -0.001576
        0.000255   -0.000143   -0.001264    0.002278    0.000000    0.000000

### `seedname_htL.dat`

INPUT. Read if $\verb#transport_mode#=\verb#lcr#$ and
$\verb#tran_read_ht#=\verb#.TRUE.#$. The file must be written in the
same way as in `seedname_htB.dat`. The first line can be any comment you
want. The second line gives `tran_num_ll`. `tran_num_ll` in
`seedname_htL.dat` must be equal to that in `seedname.win`. The code
will stop otherwise.

     Created by a WANNIER user
       105
        0.316879    0.000000   -2.762434    0.048956    0.000000   -0.016639
        0.000000    0.000000    0.000000    0.000000    0.000000   -2.809405
        .
        .
        .
        0.000000    0.078188    0.000000    0.000000   -2.086453   -0.001535
        0.007878   -0.545485  -10.525435
       105
        0.000000    0.000000    0.000315   -0.000294    0.000000    0.000085
        0.000000    0.000000    0.000000    0.000000    0.000000    0.000021
        .
        .
        .
        0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
        0.000000    0.000000    0.000000

### `seedname_htR.dat`

INPUT. Read if $\verb#transport_mode#=\verb#lcr#$ and
$\verb#tran_read_ht#=\verb#.TRUE.#$ and
$\verb#tran_use_same_lead#=\verb#.FALSE.#$. The file must be written in
the same way as in `seedname_htL.dat`. `tran_num_rr` in
`seedname_htR.dat` must be equal to that in `seedname.win`.

### `seedname_htC.dat`

INPUT. Read if $\verb#transport_mode#=\verb#lcr#$ and
$\verb#tran_read_ht#=\verb#.TRUE.#$. The first line can be any comment
you want. The second line gives `tran_num_cc`. The subsequent lines
contain `tran_num_cc`$\times$`tran_num_cc` $H_{mn}$ matrix, where the
indices $m$ and $n$ span all `tran_num_cc` WFs inside the central
conductor region. `tran_num_cc` in `seedname_htC.dat` must be equal to
that in `seedname.win`.

     Created by a WANNIER user
        99
      -10.499455   -0.541232    0.007684   -0.001624   -2.067078   -0.412188
        0.003217    0.076965    0.000522   -0.000414    0.000419   -2.122184
        .
        .
        .
       -0.003438    0.078545    0.024426    0.757343   -2.004899   -0.001632
        0.007807   -0.542983  -10.516896

### `seedname_htLC.dat`

INPUT. Read if $\verb#transport_mode#=\verb#lcr#$ and
$\verb#tran_read_ht#=\verb#.TRUE.#$. The first line can be any comment
you want. The second line gives `tran_num_ll` and `tran_num_lc` in the
given order. The subsequent lines contain
`tran_num_ll`$\times$`tran_num_lc` $H_{mn}$ matrix. The index $m$ spans
`tran_num_ll` WFs in the surface principal layer of semi-infinite left
lead which is in contact with the conductor region. The index $n$ spans
`tran_num_lc` WFs in the conductor region which have a non-negligible
interaction with the WFs in the semi-infinite left lead. Note that
`tran_num_lc` can be different from `tran_num_cc`.

     Created by a WANNIER user
       105    99
        0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
        0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
        .
        .
        .
       -0.000003    0.000009    0.000290    0.000001   -0.000007   -0.000008
        0.000053   -0.000077   -0.000069

### `seedname_htCR.dat`

INPUT. Read if $\verb#transport_mode#=\verb#lcr#$ and
$\verb#tran_read_ht#=\verb#.TRUE.#$. The first line can be any comment
you want. The second line gives `tran_num_cr` and `tran_num_rr` in the
given order. The subsequent lines contain
`tran_num_cr`$\times$`tran_num_rr` $H_{mn}$ matrix. The index $m$ spans
`tran_num_cr` WFs in the conductor region which have a non-negligible
interaction with the WFs in the semi-infinite right lead. The index $n$
spans `tran_num_rr` WFs in the surface principal layer of semi-infinite
right lead which is in contact with the conductor region. Note that
`tran_num_cr` can be different from `tran_num_cc`.

     Created by a WANNIER user
        99   105
       -0.000180    0.000023    0.000133   -0.000001    0.000194    0.000008
       -0.000879   -0.000028    0.000672   -0.000257   -0.000102   -0.000029
        .
        .
        .
        0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
        0.000000    0.000000    0.000000

### `seedname.unkg` {#sec:files_unkg}

INPUT. Read if $\verb#transport_mode#=\verb#lcr#$ and
$\verb#tran_read_ht#=\verb#.FALSE.#$. The first line is the number of
G-vectors at which the $\tilde{u}_{m\mathbf{k}}(\mathbf{G})$ are
subsequently printed. This number should always be 32 since 32 specific
$\tilde{u}_{m\mathbf{k}}$ are required. The following lines contain the
following in this order: The band index $m$, a counter on the number of
G-vectors, the integer co-efficient of the G-vector components $a,b,c$
(where $\mathbf{G}=a\mathbf{b}_1+b\mathbf{b}_2+c\mathbf{b}_3$), then the
real and imaginary parts of the corresponding
$\tilde{u}_{m\mathbf{k}}(\mathbf{G})$ at the $\Gamma$-point. We note
that the ordering in which the G-vectors and
$\tilde{u}_{m\mathbf{k}}(\mathbf{G})$ are printed is not important, but
the specific G-vectors are critical. The following example displays for
a single band, the complete set of $\tilde{u}_{m\mathbf{k}}(\mathbf{G})$
that are required. Note the G-vectors ($a,b,c$) needed.

          32
        1    1    0    0    0   0.4023306   0.0000000
        1    2    0    0    1  -0.0000325   0.0000000
        1    3    0    1    0  -0.3043665   0.0000000
        1    4    1    0    0  -0.3043665   0.0000000
        1    5    2    0    0   0.1447143   0.0000000
        1    6    1   -1    0   0.2345179   0.0000000
        1    7    1    1    0   0.2345179   0.0000000
        1    8    1    0   -1   0.0000246   0.0000000
        1    9    1    0    1   0.0000246   0.0000000
        1   10    0    2    0   0.1447143   0.0000000
        1   11    0    1   -1   0.0000246   0.0000000
        1   12    0    1    1   0.0000246   0.0000000
        1   13    0    0    2   0.0000338   0.0000000
        1   14    3    0    0  -0.0482918   0.0000000
        1   15    2   -1    0  -0.1152414   0.0000000
        1   16    2    1    0  -0.1152414   0.0000000
        1   17    2    0   -1  -0.0000117   0.0000000
        1   18    2    0    1  -0.0000117   0.0000000
        1   19    1   -2    0  -0.1152414   0.0000000
        1   20    1    2    0  -0.1152414   0.0000000
        1   21    1   -1   -1  -0.0000190   0.0000000
        1   22    1   -1    1  -0.0000190   0.0000000
        1   23    1    1   -1  -0.0000190   0.0000000
        1   24    1    1    1  -0.0000190   0.0000000
        1   25    1    0   -2  -0.0000257   0.0000000
        1   26    1    0    2  -0.0000257   0.0000000
        1   27    0    3    0  -0.0482918   0.0000000
        1   28    0    2   -1  -0.0000117   0.0000000
        1   29    0    2    1  -0.0000117   0.0000000
        1   30    0    1   -2  -0.0000257   0.0000000
        1   31    0    1    2  -0.0000257   0.0000000
        1   32    0    0    3   0.0000187   0.0000000
        2    1    0    0    0  -0.0000461   0.0000000
        .
        .
        .

### `seedname_u.mat`

OUTPUT. Written if $\verb#write_u_matrices#=\verb#.TRUE.#$. The first
line gives the date and time at which the file was created. The second
line states the number of kpoints `num_kpts` and the number of wannier
functions `num_wann` twice. The third line is empty. Then there are
`num_kpts` blocks of data, each of which starts with a line containing
the kpoint (in fractional coordinates of the reciprocal lattice vectors)
followed by `num_wann * num_wann` lines containing the matrix elements
(real and imaginary parts) of $\mathbf{U}^{(\mathbf{k})}$. The matrix
elements are in column-major order (ie, cycling over rows first and then
columns). There is an empty line between each block of data.

     written on 15Sep2016 at 16:33:46 
               64           8           8
    	 
       0.0000000000  +0.0000000000  +0.0000000000
       0.4468355787  +0.1394579978
      -0.0966033667  +0.4003934902
      -0.0007748974  +0.0011788678
      -0.0041177339  +0.0093821027
       .
       .
       .

       0.1250000000   0.0000000000  +0.0000000000
       0.4694005589  +0.0364941808
      +0.2287801742  -0.1135511138
      -0.4776782452  -0.0511719121
      +0.0142081014  +0.0006203139
       .
       .
       .

### `seedname_u_dis.mat`

OUTPUT. Written if $\verb#write_u_matrices#=\verb#.TRUE.#$ and
disentanglement is enabled. The first line gives the date and time at
which the file was created. The second line states the number of kpoints
`num_kpts`, the number of wannier functions `num_bands` and the number
of `num_bands`. The third line is empty. Then there are `num_kpts`
blocks of data, each of which starts with a line containing the kpoint
(in fractional coordinates of the reciprocal lattice vectors) followed
by `num_wann * num_bands` lines containing the matrix elements (real and
imaginary parts) of $\mathbf{U}^{\mathrm{dis}(\mathbf{k})}$. The matrix
elements are in column-major order (ie, cycling over rows first and then
columns). There is an empty line between each block of data.

     written on 15Sep2016 at 16:33:46 
               64           8          16
    	    
       0.0000000000  +0.0000000000  +0.0000000000
       1.0000000000  +0.0000000000
      +0.0000000000  +0.0000000000
      +0.0000000000  +0.0000000000
      +0.0000000000  +0.0000000000
       .
       .
       .
       
       0.1250000000   0.0000000000  +0.0000000000
       1.0000000000  +0.0000000000
      +0.0000000000  +0.0000000000
      +0.0000000000  +0.0000000000
      +0.0000000000  +0.0000000000
       .
       .
       .

## []{#chap:interpolation label="chap:interpolation"}Some notes on the interpolation

In `wannier90` v.2.1, a new flag `use_ws_distance` has been introduced
(and it is set to `.true.` by default since version v3.0). Setting it to
`.false.` reproduces the "standard" behavior of `wannier90` in v.2.0.1
and earlier, while setting it to `.true.` changes the interpolation
method as described below. In general, this allows a smoother
interpolation, helps reducing (a bit) the number of $k-$points required
for interpolation, and reproduces the band structure of large supercells
sampled at $\Gamma$ only (setting it to `.false.` produces instead flat
bands, which might instead be the intended behaviour for small molecules
carefully placed at the centre of the cell).

The core idea rests on the fact that the Wannier functions
$w_{n\bm{\mathrm{R}}}(\bm{\mathrm{r}})$ that we build from
$N\times M\times L$ $k-$points are actually periodic over a supercell of
size $N\times M\times L$, but when you use them to interpolate you want
them to be *zero* outside this supercell. In 1D it is pretty obvious
want we mean here, but in 3D what you really want that they are zero
outside the Wigner--Seitz cell of the $N\times M\times L$ superlattice.

The best way to impose this condition is to check that every real-space
distance that enters in the $R\to k$ Fourier transform is the shortest
possible among all the $N\times M\times L-$periodic equivalent copies.

If the distances were between unit cells, this would be trivial, but the
distances are between Wannier functions which are not centred on
$\bm{\mathrm{R}}=0$. Hence, when you want to consider the matrix element
of a generic operator $\bm{\mathrm{O}}$ (i.e., the Hamiltonian)
$\langle w_{i\bm{\mathrm{0}}}(\bm{\mathrm{r}})|\bm{\mathrm{O}}|w_{j\bm{\mathrm{R}}}(\bm{\mathrm{r}})\rangle$
you must take in account that the centre $\bm{\mathrm{\tau}}_i$ of
$w_{i\bm{\mathrm{0}}}(\bm{\mathrm{r}})$ may be very far away from
$\bm{\mathrm{0}}$ and the centre $\bm{\mathrm{\tau}}_j$ of
$w_{j\bm{\mathrm{R}}}(\bm{\mathrm{r}})$ may be very far away from
$\bm{\mathrm{R}}$.

There are many way to find the shortest possible distance between
$w_{i\bm{\mathrm{0}}}(\bm{\mathrm{r}})$ and
$w_{j\bm{\mathrm{R}}}(\bm{\mathrm{r}}-\bm{\mathrm{R}})$, the one used
here is to consider the distance
$\bm{\mathrm{d}}_{ij\bm{\mathrm{R}}} = \bm{\mathrm{\tau}}_i - (\bm{\mathrm{\tau}}_j+\bm{\mathrm{R}})$
and all its superlattice periodic equivalents
$\bm{\mathrm{d}}_{ij\bm{\mathrm{R}}}+ \bm{\mathrm{\tilde R}}_{nml}$,
with
$\bm{\mathrm{\tilde R}}_{nml} = (Nn\bm{\mathrm{a}}_1 + Mm\bm{\mathrm{a}}_2 + Ll\bm{\mathrm{a}}_3)$
and $n,l,m = {-L,-L+1,...0,...,L-1,L}$, with $L$ controlled by the
parameter `ws_search_size`.

Then,

1.  if
    $\bm{\mathrm{d}}_{ij\bm{\mathrm{R}}}+ \bm{\mathrm{\tilde R}}_{nml}$
    is inside the $N\times M \times L$ super-WS cell, then it is the
    shortest, take it and quit

2.  if it is outside the WS, then it is not the shortest, throw it away

3.  if it is on the border/corner of the WS then it is the shortest, but
    there are other choices of $(n,m,l)$ which are equivalent, find all
    of them

In all distance comparisons, a small but finite tolerance is considered,
which can be controlled with the parameter `ws_distance_tol`.

Because of how the Fourier transform is defined in the `wannier90` code
(not the only possible choice) it is only
$\bm{\mathrm{R}}+\bm{\mathrm{\tilde R}}_{nml}$ that enters the
exponential, but you still have to consider the distance among the
actual centres of the Wannier functions. Using the centres of the
unit-cell to which the Wannier functions belong is not enough (but is
easier, and saves you one index).

Point 3 is not stricly necessary, but using it helps enforcing the
symmetry of the system in the resulting band structure. You will get
some small but evident symmetry breaking in the band plots if you just
pick one of the equivalent $\bm{\mathrm{\tilde R}}$ vectors.

Note that in some cases, all this procedure does absolutely nothing, for
instance if all the Wannier function centres are very close to 0 (e.g.,
a molecule carefully placed in the periodic cell).

In some other cases, the effect may exist but be imperceptible. E.g., if
you use a very fine grid of $k-$points, even if you don't centre each
functions perfectly, the periodic copies will still be so far away that
the change in centre applied with $\tt use\_ws\_distance$ does not
matter.

When instead you use few $k-$points, activating the
$\tt use\_ws\_distance$ may help a lot in avoiding spurious oscillations
of the band structure even when the Wannier functions are well
converged.

## Sample Input Files {#chap:files}

### Master input file: `seedname.win` {#winfile}

    num_wann          : 4 
    mp_grid           : 4 4 4 
    num_iter          : 100
    postproc_setup    : true

    begin unit_cell_cart
    ang
    -1.61 0.00 1.61
     0.00 1.61 1.61
    -1.61 1.61 0.00
    end unit_cell_cart

    begin atoms_frac
    C   -0.125  -0.125  -0.125
    C    0.125   0.125   0.125
    end atoms_frac

    bands_plot        : true
    bands_num_points  : 100
    bands_plot_format : gnuplot

    begin kpoint_path
    L 0.50000 0.50000 0.50000 G 0.00000 0.00000 0.00000
    G 0.00000 0.00000 0.00000 X 0.50000 0.00000 0.50000
    X 0.50000 0.00000 0.50000 K 0.62500 0.25000 0.62500
    end kpoint_path

    begin projections
    C:l=0,l=1
    end projections

    begin kpoints
    0.00 0.00 0.00
    0.00 0.00 0.25
    0.00 0.50 0.50
     .
     .
     .
    0.75 0.75 0.50
    0.75 0.75 0.75
    end kpoints

### `seedname.nnkp` {#nnkp-file}

Running `wannier90` on the above input file would generate the following
`nnkp` file:

    File written on  9Feb2006 at 15:13: 9 

    calc_only_A   :  F

    begin real_lattice
      -1.612340   0.000000   1.612340
       0.000000   1.612340   1.612340
      -1.612340   1.612340   0.000000
    end real_lattice

    begin recip_lattice
      -1.951300  -1.951300   1.951300
       1.951300   1.951300   1.951300
      -1.951300   1.951300  -1.951300
    end recip_lattice

    begin kpoints
         64
      0.00000   0.00000   0.00000   
      0.00000   0.25000   0.00000   
      0.00000   0.50000   0.00000   
      0.00000   0.75000   0.00000   
      0.25000   0.00000   0.00000   
      .
      .
      .
      0.50000   0.75000   0.75000   
      0.75000   0.00000   0.75000   
      0.75000   0.25000   0.75000   
      0.75000   0.50000   0.75000   
      0.75000   0.75000   0.75000     
    end kpoints

    begin projections
       8
      -0.12500   -0.12500   -0.12500     0  1  1 
         0.000  0.000  1.000   1.000  0.000  0.000   2.00 
      -0.12500   -0.12500   -0.12500     1  1  1 
         0.000  0.000  1.000   1.000  0.000  0.000   2.00 
      -0.12500   -0.12500   -0.12500     1  2  1 
         0.000  0.000  1.000   1.000  0.000  0.000   2.00 
      -0.12500   -0.12500   -0.12500     1  3  1 
         0.000  0.000  1.000   1.000  0.000  0.000   2.00 
       0.12500    0.12500    0.12500     0  1  1 
         0.000  0.000  1.000   1.000  0.000  0.000   2.00 
       0.12500    0.12500    0.12500     1  1  1 
         0.000  0.000  1.000   1.000  0.000  0.000   2.00 
       0.12500    0.12500    0.12500     1  2  1 
         0.000  0.000  1.000   1.000  0.000  0.000   2.00 
       0.12500    0.12500    0.12500     1  3  1 
         0.000  0.000  1.000   1.000  0.000  0.000   2.00 
    end projections

    begin nnkpts
        8
      1     2      0   0   0
      1     4      0  -1   0
      1     5      0   0   0
      1    13     -1   0   0
      1    17      0   0   0
      1    22      0   0   0
      1    49      0   0  -1
      1    64     -1  -1  -1
      2     1      0   0   0
      2     3      0   0   0
      2     6      0   0   0
      2    14     -1   0   0
      2    18      0   0   0
      2    23      0   0   0
      2    50      0   0  -1
      2    61     -1   0  -1
      .
      .
      .
     64     1      1   1   1
     64    16      0   0   1
     64    43      0   0   0
     64    48      0   0   0
     64    52      1   0   0
     64    60      0   0   0
     64    61      0   1   0
     64    63      0   0   0
    end nnkpts

    begin exclude_bands 
       4 
       1 
       2 
       3
       4
    end exclude_bands

# `postw90.x`

## Overview of the `berry` module {#ch:berry}

The `berry` module of `postw90` is called by setting ` berry = true` and
choosing one or more of the available options for `berry_task`. The
routines in the `berry` module which compute the $k$-space Berry
curvature, orbital magnetization and spin Hall conductivity are also
called when `kpath = true` and `kpath_task = {curv,morb,shc}`, or when
`kslice = true` and `kslice_task = {curv,morb,shc}`.

### Background: Berry connection and curvature

The Berry connection is defined in terms of the cell-periodic Bloch
states $\vert u_{n{\bf k}}\rangle=e^{-i{\bf k}\cdot{\bf r}}\vert
\psi_{n{\bf k}}\rangle$ as
$${\bf A}_n({\bf k})=\langle u_{n{\bf k}}\vert i\bm{\nabla}_{\bf k}\vert
u_{n{\bf k}}\rangle,$$ and the Berry curvature is the curl of the
connection,
$$\bm{\Omega}_n({\bf k})=\bm{\nabla}_{\bf k}\times {\bf A}_n({\bf k})=
-{\rm Im}
\langle \bm{\nabla}_{\bf k} u_{n{\bf k}}\vert \times
\vert\bm{\nabla}_{\bf k} u_{n{\bf k}}\rangle.$$ These two quantities
play a central role in the description of several electronic properties
of crystals [@xiao-rmp10]. In the following we will work with a matrix
generalization of the Berry connection,
$${\bf A}_{nm}({\bf k})=\langle u_{n{\bf k}}\vert i\bm{\nabla}_{\bf k}\vert
u_{m{\bf k}}\rangle={\bf A}_{mn}^*({\bf k}),
\label{eq:berry-connection-matrix}$$ and write the curvature as an
antisymmetric tensor, $$\label{eq:curv}
\Omega_{n,\alpha\beta}({\bf k}) =\epsilon_{\alpha\beta\gamma}
\Omega_{n,\gamma}({\bf k})=-2{\rm Im}\langle 
\nabla_{k_\alpha} u_{n\bf k}\vert \nabla_{k_\beta} u_{n\bf k}\rangle.$$

### `berry_task=kubo`: optical conductivity and joint density of states 

The Kubo-Greenwood formula for the optical conductivity of a crystal in
the independent-particle approximation reads
$$\sigma_{\alpha\beta}(\hbar\omega)=\frac{ie^2\hbar}{N_k\Omega_c}
\sum_{\bf k}\sum_{n,m}
\frac{f_{m{\bf k}}-f_{n{\bf k}}}
     {\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}}
\frac{\langle\psi_{n{\bf k}}\vert v_\alpha\vert\psi_{m{\bf k}}\rangle
      \langle\psi_{m{\bf k}}\vert v_\beta\vert\psi_{n{\bf k}}\rangle}
{\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}-(\hbar\omega+i\eta)}.$$
Indices $\alpha,\beta$ denote Cartesian directions, $\Omega_c$ is the
cell volume, $N_k$ is the number of $k$-points used for sampling the
Brillouin zone, and $f_{n{\bf k}}=f(\varepsilon_{n{\bf k}})$ is the
Fermi-Dirac distribution function. $\hbar\omega$ is the optical
frequency, and $\eta>0$ is an adjustable smearing parameter with units
of energy.

The off-diagonal velocity matrix elements can be expressed in terms of
the connection matrix [@blount-ssp62], $$\label{eq:velocity_mat}
\langle\psi_{n{\bf k}}\vert {\bf v} \vert\psi_{m{\bf k}}\rangle=
-\frac{i}{\hbar}(\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}})
{\bf A}_{nm}({\bf k})\,\,\,\,\,\,\,\,(m\not= n).$$ The conductivity
becomes $$\begin{aligned}
\label{eq:sig-bz}
\sigma_{\alpha\beta}(\hbar\omega)&=
\frac{1}{N_k}\sum_{\bf k}\sigma_{{\bf k},\alpha\beta}(\hbar\omega)\\
\label{eq:sig-k}
\sigma_{{\bf k},\alpha\beta}(\hbar\omega)&=\frac{ie^2}{\hbar\Omega_c}\sum_{n,m}
(f_{m{\bf k}}-f_{n{\bf k}})
\frac{\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}}
{\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}-(\hbar\omega+i\eta)}
A_{nm,\alpha}({\bf k})A_{mn,\beta}({\bf k}).
\end{aligned}$$

Let us decompose it into Hermitian (dissipative) and anti-Hermitean
(reactive) parts. Note that $$\label{eq:lorentzian}
\overline{\delta}(\varepsilon)=\frac{1}{\pi}{\rm Im}
\left[\frac{1}{\varepsilon-i\eta}\right],$$ where $\overline{\delta}$
denotes a "broadended" delta-function. Using this identity we find for
the Hermitean part $$\label{eq:sig-H}
\sigma_{{\bf k},\alpha\beta}^{\rm H}(\hbar\omega)=-\frac{\pi e^2}{\hbar\Omega_c}
\sum_{n,m}(f_{m{\bf k}}-f_{n{\bf k}})
(\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}})
A_{nm,\alpha}({\bf k})A_{mn,\beta}({\bf k})
\overline{\delta}(\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}-\hbar\omega).$$
Improved numerical accuracy can be achieved by replacing the Lorentzian
([\[eq:lorentzian\]](#eq:lorentzian){reference-type="ref"
reference="eq:lorentzian"}) with a Gaussian, or other shapes. The
analytical form of $\overline{\delta}(\varepsilon)$ is controlled by the
keyword `[kubo_]smr_type`.

The anti-Hermitean part of
Eq. ([\[eq:sig-k\]](#eq:sig-k){reference-type="ref"
reference="eq:sig-k"}) is given by $$\label{eq:sig-AH}
\sigma_{{\bf k},\alpha\beta}^{\rm AH}(\hbar\omega)=\frac{ie^2}{\hbar\Omega_c}
\sum_{n,m}(f_{m{\bf k}}-f_{n{\bf k}})
{\rm Re}\left[ \frac{\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}}
                    {\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}
                     -(\hbar\omega+i\eta)}
\right]
A_{nm,\alpha}({\bf k})A_{mn,\beta}({\bf k}).$$ Finally the joint density
of states is $$\label{eq:jdos}
\rho_{cv}(\hbar\omega)=\frac{1}{N_k}\sum_{\bf k}\sum_{n,m}
f_{n{\bf k}}(1-f_{m{\bf k}})
\overline{\delta}(\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}}-\hbar\omega).$$

Equations ([\[eq:lorentzian\]](#eq:lorentzian){reference-type="ref"
reference="eq:lorentzian"}--[\[eq:jdos\]](#eq:jdos){reference-type="ref"
reference="eq:jdos"}) contain the parameter $\eta$, whose value can be
chosen using the keyword\
` [kubo_]smr_fixed_en_width`. Better results can often be achieved by
adjusting the value of $\eta$ for each pair of states, i.e.,
$\eta\rightarrow \eta_{nm\bf k}$. This is done as follows (see
description of the keyword `adpt_smr_fac`)
$$\eta_{nm{\bf k}}=\alpha\vert \bm{\nabla}_{\bf k}
(\varepsilon_{m{\bf k}}-\varepsilon_{n{\bf k}})\vert \Delta k.$$

The energy eigenvalues $\varepsilon_{n\bf k}$, band velocities
$\bm{\nabla}_{\bf k}\varepsilon_{n{\bf k}}$, and off-diagonal Berry
connection ${\bf A}_{nm}({\bf k})$ entering the previous four equations
are evaluated over a $k$-point grid by Wannier interpolation, as
described in Refs. [@wang-prb06; @yates-prb07]. After averaging over the
Brillouin zone, the Hermitean and anti-Hermitean parts of the
conductivity are assembled into the symmetric and antisymmetric tensors
$$\begin{aligned}
\sigma^{\rm S}_{\alpha\beta}&=
{\rm Re}\sigma^{\rm H}_{\alpha\beta}+i{\rm Im}\sigma^{\rm AH}_{\alpha\beta}\\
\sigma^{\rm A}_{\alpha\beta}&=
{\rm Re}\sigma^{\rm AH}_{\alpha\beta}+i{\rm Im}\sigma^{\rm H}_{\alpha\beta},
\end{aligned}$$ whose independent components are written as a function
of frequency onto nine separate files.

### `berry_task=ahc`: anomalous Hall conductivity

The antisymmetric tensor $\sigma^{\rm A}_{\alpha\beta}$ is odd under
time reversal, and therefore vanishes in non-magnetic systems, while in
ferromagnets with spin-orbit coupling it is generally nonzero. The
imaginary part ${\rm Im}\sigma^{\rm H}_{\alpha\beta}$ describes magnetic
circular dichroism, and vanishes as $\omega\rightarrow
0$. The real part ${\rm Re}\sigma^{\rm AH}_{\alpha\beta}$ describes the
anomalous Hall conductivity (AHC), and remains finite in the static
limit.

The intrinsic dc AHC is obtained by setting $\eta=0$ and $\omega=0$ in
Eq. ([\[eq:sig-AH\]](#eq:sig-AH){reference-type="ref"
reference="eq:sig-AH"}). The contribution from point ${\bf k}$ in the
Brillouin zone is
$$\sigma^{\rm AH}_{{\bf k},\alpha\beta}(0)=\frac{2e^2}{\hbar\Omega_c}
\sum_{n,m}f_{n\bf k}(1-f_{m\bf k})
{\rm Im}\langle \nabla_{k_\alpha} u_{n\bf k}\vert u_{m\bf k}\rangle
\langle u_{m\bf k}\vert\nabla_{k_\beta} u_{n\bf k}\rangle,$$ where we
replaced $f_{n\bf k}-f_{m\bf k}$ with
$f_{n\bf k}(1-f_{m\bf k})-f_{m\bf k}(1-f_{n\bf k})$.

This expression is not the most convenient for *ab initio* calculations,
as the sums run over the complete set of occupied and empty states. In
practice the sum over empty states can be truncated, but a relatively
large number should be retained to obtain accurate results. Using the
resolution of the identity $1=\sum_m \vert u_{m\bf
  k}\rangle \langle u_{m\bf k}\vert$ and noting that the term
$\sum_{n,m}f_{n\bf k}f_{m\bf k}(\ldots)$ vanishes identically, we arrive
at the celebrated formula for the intrinsic AHC in terms of the Berry
curvature, $$\begin{aligned}
\label{eq:ahc}
\sigma^{\rm AH}_{\alpha\beta}(0)&=\frac{e^2}{\hbar}
\frac{1}{N_k\Omega_c}\sum_{\bf k}(-1)\Omega_{\alpha\beta}({\bf k}),\\
%\sum_n (-1)f_{n\bf k}\Omega_{n,\alpha\beta}({\bf k}).
\label{eq:curv-occ}
\Omega_{\alpha\beta}({\bf k})&=\sum_n f_{n\bf k}\Omega_{n,\alpha\beta}({\bf k}).
\end{aligned}$$ Note that only *occupied* states enter this expression.
Once we have a set of Wannier functions spanning the valence bands
(together with a few low-lying conduction bands, typically)
Eq. ([\[eq:ahc\]](#eq:ahc){reference-type="ref" reference="eq:ahc"}) can
be evaluated by Wannier interpolation as described in
Refs. [@wang-prb06; @lopez-prb12], with no truncation involved.

### `berry_task=morb`: orbital magnetization

The ground-state orbital magnetization of a crystal is given
by [@xiao-rmp10; @ceresoli-prb06] $$\begin{aligned}
\label{eq:morb}
{\bf M}^{\rm orb}&=\frac{e}{\hbar}
%\int_{\rm BZ}\frac{d{\bf k}}{(2\pi)^3}
\frac{1}{N_k\Omega_c}\sum_{\bf k}{\bf M}^{\rm orb}({\bf k}),\\
\label{eq:morb-k}
{\bf M}^{\rm orb}({\bf k})&=
\sum_n\,\frac{1}{2}f_{n{\bf k}}\,
{\rm Im}\,\langle \bm{\nabla}_{\bf k}u_{n{\bf k}}\vert
\times
\left(H_{\bf k}+\varepsilon_{n{\bf k}}-2\varepsilon_F\right)
\vert \bm{\nabla}_{\bf k}u_{n{\bf k}}\rangle,
\end{aligned}$$ where $\varepsilon_F$ is the Fermi energy. The
Wannier-interpolation calculation is described in Ref. [@lopez-prb12].
Note that the definition of ${\bf M}^{\rm orb}({\bf k})$ used here
differs by a factor of $-1/2$ from the one in Eq. (97) and Fig. 2 of
that work.

### `berry_task=shc`: spin Hall conductivity {#sec:shc}

The Kubo-Greenwood formula for the intrinsic spin Hall conductivity
(SHC) of a crystal in the independent-particle approximation reads
[@qiao-prb2018; @ryoo-prb2019; @guo-prl2008] $$\label{eq:kubo_shc}
\sigma_{\alpha\beta}^{\text{spin}\gamma}(\omega) =  \frac{\hbar}{\Omega_c N_k}
\sum_{\bm{k}}\sum_{n} f_{n\bm{k}} \\
\sum_{m \neq n}
\frac{2\operatorname{Im}[\langle n\bm{k}| \hat{j}_{\alpha}^{\gamma}|m\bm{k}\rangle
	\langle m\bm{k}| -e\hat{v}_{\beta}|n\bm{k}\rangle]}
{(\epsilon_{n\bm{k}}-\epsilon_{m\bm{k}})^2-(\hbar\omega +i\eta)^2}.$$
The spin current operator $\hat{j}_{\alpha}^{\gamma}=
\frac{1}{2}\{\hat{s}_{\gamma},\hat{v}_{\alpha}\}$ where the spin
operator $\hat{s}_{\gamma}=\frac{\hbar}{2}\hat{\sigma}_{\gamma}$.
Indices $\alpha,\beta$ denote Cartesian directions, $\gamma$ denotes the
direction of spin, commonly $\alpha = x, \beta = y, \gamma = z$.
$\Omega_c$ is the cell volume, $N_k$ is the number of $k$-points used
for sampling the Brillouin zone, and
$f_{n{\bf k}}=f(\varepsilon_{n{\bf k}})$ is the Fermi-Dirac distribution
function. $\hbar\omega$ is the optical frequency, and $\eta>0$ is an
adjustable smearing parameter with unit of energy.

The velocity matrix element in the numerator is the same as
Eq. ([\[eq:velocity_mat\]](#eq:velocity_mat){reference-type="ref"
reference="eq:velocity_mat"}), so the only unknown quantity is the spin
current matrix
$\langle n\bm{k}| \hat{j}_{\alpha}^{\gamma}|m\bm{k}\rangle$. We can use
Wannier interpolation technique to efficiently calculate this matrix,
and there are two derivation according to the degree of approximation. A
noteworthy difference is the way in which two *ab-initio* matrix
elements are evaluated,
$$\langle u_{n{\bf k}}\vert\sigma_\gamma H_{\bf k}\vert u_{m{\bf k}+{\bf b}}\rangle, \langle u_{n{\bf k}}\vert\sigma_\gamma \vert u_{m{\bf k}+{\bf b}}\rangle, \gamma = x, y, z$$
These are evaluated by `pw2wannier90` using Ryoo's method. In contrast,
Qiao's method does not require `pw2wannier90`, but it assumes an
approximation
$1\approx\sum_{ l\in ab-initio{\rm \,bands}}|u_{l\bm{k}}\rangle \langle u_{l\bm{k}}|$.
You can choose which method to evaluate this value with `shc_method` in
the input file. For a full derivation please refer to
Ref. [@qiao-prb2018] or Ref. [@ryoo-prb2019].

The Eq. ([\[eq:kubo_shc\]](#eq:kubo_shc){reference-type="ref"
reference="eq:kubo_shc"}) can be further separated into band-projected
Berry curvature-like term $$\label{eq:kubo_shc_berry}
\Omega_{n,\alpha\beta}^{\text{spin}\gamma}(\bm{k}) = {\hbar}^2 \sum_{
	m\ne n}\frac{-2\operatorname{Im}[\langle n\bm{k}| 
	\frac{1}{2}\{\hat{\sigma}_{\gamma},\hat{v}_{\alpha}\}|m\bm{k}\rangle
	\langle m\bm{k}| \hat{v}_{\beta}|n\bm{k}\rangle]}
{(\epsilon_{n\bm{k}}-\epsilon_{m\bm{k}})^2-(\hbar\omega+i\eta)^2},$$
$k$-resolved term which sums over occupied bands
$$\label{eq:kubo_shc_berry_sum}
\Omega_{\alpha\beta}^{\text{spin}\gamma}(\bm{k}) = \sum_{n}
f_{n\bm{k}} \Omega_{n,\alpha\beta}^{\text{spin}\gamma}(\bm{k}),$$ and
the SHC is $$\sigma_{\alpha\beta}^{\text{spin}\gamma}(\omega) = 
-\frac{e^2}{\hbar}\frac{1}{\Omega_c N_k}\sum_{\bm{k}}
\Omega_{\alpha\beta}^{\text{spin}\gamma}(\bm{k}).$$ The unit of the
$\Omega_{n,\alpha\beta}^{\text{spin}\gamma}(\bm{k})$ is
$\text{length}^{2}$ (Angstrom$^2$ or Bohr$^2$, depending on your choice
of `berry_curv_unit` in the input file), and the unit of
$\sigma_{\alpha\beta}^{\text{spin}\gamma}$ is $(\hbar/e)$S/cm (the unit
is written in the header of the output file). The case of $\omega=0$
corresponds to direct current (dc) SHC while that of $\omega\ne0$
corresponds to alternating current (ac) SHC or frequency-dependent SHC.
Note in some papers
Eq. ([\[eq:kubo_shc_berry\]](#eq:kubo_shc_berry){reference-type="ref"
reference="eq:kubo_shc_berry"}) is called as spin Berry curvature.
However, it was pointed out by Ref. [@Gradhand_2012] that this name is
misleading, so we use a somewhat awkward name "Berry curvature-like
term" to refer to
Eq. ([\[eq:kubo_shc_berry\]](#eq:kubo_shc_berry){reference-type="ref"
reference="eq:kubo_shc_berry"}). The $k$-resolved term
Eq. ([\[eq:kubo_shc_berry_sum\]](#eq:kubo_shc_berry_sum){reference-type="ref"
reference="eq:kubo_shc_berry_sum"}) can be used to draw `kslice` plot,
and the band-projected Berry curvature-like term
Eq. ([\[eq:kubo_shc_berry\]](#eq:kubo_shc_berry){reference-type="ref"
reference="eq:kubo_shc_berry"}) can be used to color the `kpath` plot.

Same as the case of optical conductivity, the parameter $\eta$ contained
in the
Eq. ([\[eq:kubo_shc_berry\]](#eq:kubo_shc_berry){reference-type="ref"
reference="eq:kubo_shc_berry"}) can be chosen using the keyword
`[kubo_]smr_fixed_en_width`. Also, adaptive smearing can be employed by
the keyword `[kubo_]adpt_smr` (see Examples 29 and 30 in the Tutorial).

Please cite the following paper [@qiao-prb2018] or  [@ryoo-prb2019] when
publishing SHC results obtained using this method:

> Junfeng Qiao, Jiaqi Zhou, Zhe Yuan, and Weisheng Zhao,\
> *Calculation of intrinsic spin Hall conductivity by Wannier
> interpolation*,\
> Phys. Rev. B. 98, 214402 (2018), DOI:10.1103/PhysRevB.98.214402.

or

> Ji Hoon Ryoo, Cheol-hwan Park, and Ivo Souza,\
> *Computation of intrinsic spin Hall conductivities from first
> principles using maximally localized Wannier functions*,\
> Phys. Rev. B. 99, 235113 (2019), DOI:10.1103/PhysRevB.99.235113.

### `berry_task=sc`: shift current

The shift-current contribution to the second-order response is
characterized by a frequency-dependent third-rank tensor [@sipe-prb00]
$$\label{eq:shiftcurrent}
\begin{split}
\sigma^{abc}(0;\omega,-\omega)=&-\frac{i\pi e^3}{4\hbar^2 \Omega_c N_k}
\sum_{\bm{k}} \sum_{n,m}(f_{n\bm{k}}-f_{m\bm{k}})
\times
\left(r^b_{ mn}(\bm{k})r^{c;a}_{nm}(\bm{k}) + r^c_{mn}(\bm{k})r^{b;a}_{ nm}(\bm{k})\right)\\
&\times \left[\delta(\omega_{mn\bm{k}}-\omega)+\delta(\omega_{nm\bm{k}}-\omega)\right],
\end{split}$$ where $a,b,c$ are spatial indexes and
$\omega_{mn\bm{k}}=(\epsilon_{n\bm{k}}-\epsilon_{m\bm{k}})/\hbar$. The
expression in
Eq. [\[eq:shiftcurrent\]](#eq:shiftcurrent){reference-type="ref"
reference="eq:shiftcurrent"} involves the dipole matrix element
$$\label{eq:r}
r^a_{ nm}(\bm{k})=(1-\delta_{nm})A^a_{ nm}(\bm{k}),$$ and its
*generalized derivative* $$\label{eq:gen-der}
r^{a;b}_{nm}(\bm{k})=\partial_{k_{b}} r^a_{nm}(\bm{k})
-i\left(A^b_{nn}(\bm{k})-A^b_{ mm}(\bm{k})\right)r^a_{ nm}(\bm{k}).$$
The first-principles evaluation of the above expression is technically
challenging due to the presence of an extra $k$-space derivative. The
implementation in `wannier90` follows the scheme proposed in
Ref. [@ibanez-azpiroz_ab_2018], following the spirit of the
Wannier-interpolation method for calculating the AHC [@wang-prb06] by
reformulating $k\cdot p$ perturbation theory within the subspace of
wannierized bands. This strategy inherits the practical advantages of
the sum-over-states approach, but without introducing the truncation
errors usually associated with this procedure [@sipe-prb00].

As in the case of the optical conductivity, a broadened delta function
can be applied in
Eq. [\[eq:shiftcurrent\]](#eq:shiftcurrent){reference-type="ref"
reference="eq:shiftcurrent"} by means of the parameter $\eta$ (see
Eq. [\[eq:lorentzian\]](#eq:lorentzian){reference-type="ref"
reference="eq:lorentzian"}) using the keyword
`[kubo_]smr_fixed_en_width`, and adaptive smearing can be employed using
the keyword `[kubo_]adpt_smr`.

Please cite Ref. [@ibanez-azpiroz_ab_2018] when publishing shift-current
results using this method.

### `berry_task=kdotp`: $k\cdot p$ coefficients {#sec:kdotp}

Consider a Hamiltonian $$\label{eq:H}
H=H^{0}+H^{\prime}$$ where the eigenvalues $E_{n}$ and eigenfunctions
$\vert n\rangle$ of $H^{0}$ are known, and $H^{\prime}$ is a
perturbation. In a nutshell, quasi-degenerate perturbation theory
assumes that the set of eigenfunctions of $H^0$ can be divided into
subsets $A$ and $B$ that are weakly coupled by $H^{\prime}$, and that we
are only interested in subset $A$. This theory asserts that a
transformed Hamiltonian $\tilde{H}$ exists within subspace $A$ such that
$$\label{eq:pert-exp}
\tilde{H}=\tilde{H}^{0}+\tilde{H}^{1}+\tilde{H}^{2} + \cdots$$ where
$\tilde{H}^{j}$ contain matrix elements of $H^{\prime}$ to the $j$th
power. According to Appendix B of Ref [@winkler_spin-orbit_2003], the
first three terms are $$\begin{aligned}
\label{eq:pert-matelem0}
& \tilde{H}^{0}_{mm'} = H^{0}_{mm'},\\
\label{eq:pert-matelem1}
& \tilde{H}^{1}_{mm'} = H^{'}_{mm'},\\
\label{eq:pert-matelem2}
& \tilde{H}^{2}_{mm'} = \dfrac{1}{2}\sum_{l}H^{'}_{ml}H^{'}_{lm'}
\left( 
\dfrac{1}{E_{m}-E_{l}}+\dfrac{1}{E_{m'}-E_{l}}
\right),
\end{aligned}$$ where $m,m'\in A$ and $l\in B$. The approximation
$\tilde{H}\sim \tilde{H}^{0}+\tilde{H}^{1}$ amounts to truncating $H$ to
the $A$ subspace. By further including $\tilde{H}^{2}$, the coupling to
the $B$ subspace is incorporated approximately, "renormalizing" the
elements of the truncated matrix.

We adopt the notation described in Sec. III.B of Ref. [@wang-prb06]. We
shift the origin of $k$ space to the point where the band edge (or some
other band extremum of interest) is located, and Taylor expand around
that point the Wannier-gauge Hamiltonian, $$\label{eq:HW-exp}
H^{(W)}(\bm{k})=H^{(W)}(0)
+\sum_{a}H_{a}^{(W)}(0)k_{a}
+\dfrac{1}{2}\sum_{ab}H_{ab}^{(W)}(0)k_{a}k_{b}
+ \mathcal{O}(k^{3})$$ where $a,b=x,y,z$, and $$\begin{aligned}
&H_{a}^{(W)}(0)=\left. \dfrac{\partial H^{(W)}(\bm{k})}{\partial k_{a}}\right\rvert_{\bm{k}=0}\\
&H_{ab}^{(W)}(0)=\left. \dfrac{\partial^{2} H^{(W)}(\bm{k})}{\partial k_{a}\partial k_{b}}\right\rvert_{\bm{k}=0}
\end{aligned}$$

We now apply to $H^{(W)}(\bm{k})$ a similarity transformation $U(0)$
that diagonalizes $H^{(W)}(0)$, and call the transformed Hamiltonian
$H(\bm{k})$, $$\label{eq:Hbar}
H(\bm{k})=\overbrace{\overline{H}}^{H^{0}} + \overbrace{\sum_{a}\overline{H}_{a}k_{a}
+\dfrac{1}{2}\sum_{ab}\overline{H}_{ab}k_{a}k_{b}}^{H^{\prime}} + \mathcal{O}(k^{3}),$$
where we introduced the notation
$$\overline{\mathcal{O}}=U^{\dagger}(0)\mathcal{O}^{(W)}(0)U(0),$$ and
applied it to $\mathcal{O}=H,{H}_{a},{H}_{ab}$. We can now apply
quasi-degenerate perturbation theory by choosing the diagonal matrix
$\overline{H}$ as our $H^{0}$, and the remaining (nondiagonal) terms in
Eq. [\[eq:Hbar\]](#eq:Hbar){reference-type="ref" reference="eq:Hbar"} as
$H^{\prime}$. Collecting terms in
Eq. ([\[eq:pert-exp\]](#eq:pert-exp){reference-type="ref"
reference="eq:pert-exp"}) up to second order in $k$ we get
$$\label{eq:Htilde}
 \tilde{H}_{mm'}(\bm{k}) = 
\overline{H}_{mm'} + \sum_{a} \left(\overline{H}_{a}\right)_{mm'}k_{a}
 + \dfrac{1}{2}\sum_{a,b}\left[
\left(\overline{H}_{ab}\right)_{mm'} + \left({T}_{ab}\right)_{mm'}
\right]k_{a}k_{b}+ \mathcal{O}(k^{3}),$$ where $m,m'\in A$ and we have
defined the virtual-transition matrix $$\label{eq:Tab}
\left({T}_{ab}\right)_{mm'}=\sum_{l\in B}
\left(\overline{H}_{a}\right)_{ml}\left(\overline{H}_{b}\right)_{lm'} 
 \times
\left( 
\dfrac{1}{E_{m}-E_{l}}+\dfrac{1}{E_{m'}-E_{l}}
\right)
= 
\left({T}_{ab}\right)_{m'm}^{*}.$$ (The $T_{ab}$ term in Eq.
[\[eq:Htilde\]](#eq:Htilde){reference-type="ref" reference="eq:Htilde"}
gives an Hermitean contribution to $\tilde{H}(\bm{k})$ only after
summing over $a$ and $b$, whereas the other terms are Hermitean already
before summing.)

The implementation in `wannier90` follows the scheme proposed in
Ref. [@ibanez-azpiroz-ArXiv2019], and outputs $\overline{H}_{mm'}$ in
`seedname-kdotp_0.dat`, $\left(\overline{H}_{a}\right)_{mm'}$ in
`seedname-kdotp_1.dat`, and
$\left[\left(\overline{H}_{ab}\right)_{mm'} + \left({T}_{ab}\right)_{mm'}\right]/2$
in `seedname-kdotp_2.dat`.

Please cite Ref. [@ibanez-azpiroz-ArXiv2019] when publishing $k\cdot p$
results using this method.

### Needed matrix elements

All the quantities entering the formulas for the optical conductivity
and AHC can be calculated by Wannier interpolation once the Hamiltonian
and position matrix elements $\langle {\bf 0}n\vert H\vert
{\bf R}m\rangle$ and $\langle {\bf 0}n\vert {\bf r}\vert {\bf
  R}m\rangle$ are known [@wang-prb06; @yates-prb07]. Those matrix
elements are readily available at the end of a standard MLWF calculation
with `wannier90`. In particular, $\langle {\bf
  0}n\vert {\bf r}\vert {\bf R}m\rangle$ can be calculated by Fourier
transforming the overlap matrices in Eq. (1.7),
$$\langle u_{n{\bf k}}\vert u_{m{\bf k}+{\bf b}}\rangle.$$ Further
Wannier matrix elements are needed for the orbital
magnetization [@lopez-prb12]. In order to calculate them using Fourier
transforms, one more piece of information must be taken from the
$k$-space *ab-initio* calculation, namely, the matrices
$$\langle u_{n{\bf k}+{\bf b}_1}\vert
H_{\bf k}\vert u_{m{\bf k}+{\bf b}_2}\rangle$$ over the *ab-initio*
$k$-point mesh [@lopez-prb12]. These are evaluated by `pw2wannier90`,
the interface routine between ` pwscf` and `wannier90`, by adding to the
input file ` seedname.pw2wan` the line $${\tt
%\begin{quote}
write\_uHu = .true.
%\end{quote}
}$$ The calculation of spin Hall conductivity needs the spin matrix
elements
$$\langle u_{n{\bf k}}\vert \sigma_\gamma \vert u_{m{\bf k}}\rangle, 
\gamma = x, y, z$$ from the *ab-initio* $k$-point mesh. These are also
evaluated by `pw2wannier90` by adding to the input file
` seedname.pw2wan` the line $${\tt
	%\begin{quote}
	write\_spn = .true.
	%\end{quote}
}$$ If one uses Ryoo's method to calculate spin Hall conductivity, the
further matrix elements are needed: $$\langle u_{n{\bf k}}\vert
\sigma_\gamma H_{\bf k}\vert u_{m{\bf k}+{\bf b}}\rangle, \langle u_{n{\bf k}}\vert
\sigma_\gamma \vert u_{m{\bf k}+{\bf b}}\rangle,
\gamma = x, y, z$$ and these are evaluated by adding to the input file
` seedname.pw2wan` the lines $${\tt
	write\_sHu = .true.
}$$ $${\tt
	write\_sIu = .true.
}$$

## Overview of the `gyrotropic` module []{#ch:gyrotropic label="ch:gyrotropic"}

The `gyrotropic` module of `postw90` is called by setting
` gyrotropic = true` and choosing one or more of the available options
for `gyrotropic_task`. The module computes the quantities, studied in
[@tsirkin-arxiv17], where more details may be found.

### `gyrotropic_task=-d0`: the Berry curvature dipole 

The traceless dimensionless tensor $$\label{eq:D_ab}
D_{ab}=\int[d{\bm k}]\sum_n
\frac{\partial E_n}{\partial{k_a}}
\Omega_n^b
\left(-\frac{\partial f_0}{\partial E}\right)_{E=E_n},$$

### `gyrotropic_task=-dw`: the finite-frequency generalization of the Berry curvature dipole 

$$\label{eq:D-tilde}
\widetilde{D}_{ab}(\omega)=\int[d{\bm k}]\sum_n
\frac{\partial E_n}{\partial{k_a}}\widetilde\Omega^b_n(\omega)
\left(-\frac{\partial f_0}{\partial E}\right)_{E=E_n},$$

where $\widetilde{\bm\Omega}_{{\bm k}n}(\omega)$ is a finite-frequency
generalization of the Berry curvature: $$\label{eq:curv-w}
\widetilde{\bm\Omega}_{{\bm k}n}(\omega)=-
\sum_m\,\frac{\omega^2_{{\bm k}mn}}{\omega^2_{{\bm k}mn}-\omega^2}
{\rm Im}\left({\bm A}_{{\bm k}nm}\times{\bm A}_{{\bm k}mn}\right)$$
Contrary to the Berry curvature, the divergence of
$\tilde{\bm\Omega}_{{\bm k}n}(\omega)$ is generally nonzero. As a
result, $\widetilde{D}(\omega)$ can have a nonzero trace at finite
frequencies, $\tilde{D}_\|\neq-2\tilde{D}_\perp$ in Te.

### `gyrotropic_task=-C`: the ohmic conductivity 

In the constant relaxation-time approximation the ohmic conductivity is
expressed as $\sigma_{ab}=(2\pi e\tau/\hbar)C_{ab}$, with
$$\label{eq:C_ab}
C_{ab}=\frac{e}{h}\int[d{\bm k}]\sum_n\,
\frac{\partial E_n}{\partial{k_a}} \frac{\partial E_n}{\partial{k_b}}
\left(-\frac{\partial f_0}{\partial E}\right)_{E=E_n}$$ a positive
quantity with units of surface current density (A/cm).

### `gyrotropic_task=-K`: the kinetic magnetoelectric effect (kME) 

A microscopic theory of the intrinsic kME effect in bulk crystals was
recently developed [@yoda-sr15; @zhong-prl16].

The response is described by $$\label{eq:K_ab}
K_{ab}=\int[d{\bm k}]\sum_n\frac{\partial E_n}{\partial{k_a}} m_n^b 
\left(-\frac{\partial f_0}{\partial E}\right)_{E=E_n},$$ which has the
same form as Eq. ([\[eq:D_ab\]](#eq:D_ab){reference-type="ref"
reference="eq:D_ab"}), but with the Berry curvature replaced by the
intrinsic magnetic moment ${\bm m}_{{\bm k}n}$ of the Bloch electrons,
which has the spin and orbital components given by [@xiao-rmp10]
$$\begin{aligned}
\label{eq:m-spin}
m^{\rm spin}_{{\bm k}n}&=&-\frac{1}{2}g_s\mu_{\rm B} \langle\psi_{{\bm k}
      n}\vert\bf \sigma\vert\psi_{{\bm k}n}\rangle\\
\label{eq:m-orb}
{\bm m}^{\rm orb}_{{\bm k}n}&=&\frac{e}{2\hbar}{\rm Im}
\langle{\bm\partial}_{\bm k}u_{{\bm k}n}\vert\times
(H_{\bm k}-E_{{\bm k}n})\vert{\bm\partial}_{\bm k}u_{{\bm k}n}\rangle,
\end{aligned}$$ where $g_s\approx 2$ and we chose $e>0$.

### `gyrotropic_task=-dos`: the density of states 

The density of states is calculated with the same width and type of
smearing, as the other properties of the `gyrotropic` module

### `gyrotropic_task=-noa`: the interband contributionto the natural optical activity 

Natural optical rotatory power is given by [@ivchenko-spss75]
$$\label{eq:rho-c}
\rho_0(\omega)=\frac{\omega^2}{2c^2}{\rm Re}\,\gamma_{xyz}(\omega).$$
for light propagating ling the main symmetry axis of a crystal $z$. Here
$\gamma_{xyz}(\omega)$ is an anti-symmetric (in $xy$) tensor with units
of length, which has both inter- and intraband contributions.

Following Ref. [@malashevich-prb10] for the interband contribution we
writewe write, with $\partial_c\equiv\partial/\partial k_c$,
$$\begin{gathered}
{\rm Re}\,\gamma_{abc}^{\mathrm{inter}}(\omega)=\frac{e^2}{\varepsilon_0\hbar^2}
\int[d{\bm k}]
\sum_{n,l}^{o,e}\,
\Bigl[ \frac{1}{\omega_{ln}^2-\omega^2} 
{\rm Re}\left(A_{ln}^bB_{nl}^{ac}-A_{ln}^aB_{nl}^{bc}\right) \\
-\frac{3\omega_{ln}^2-\omega^2}{(\omega_{ln}^2-\omega^2)^2} 
\partial_c(E_l+E_n){\rm Im}\left(A_{nl}^aA_{ln}^b\right)   
\Bigr].
\label{eq:gamma-inter}
\end{gathered}$$ The summations over $n$ and $l$ span the occupied ($o$)
and empty ($e$) bands respectively, $\omega_{ln}=(E_l-E_n)/\hbar$, and
${\bm A}_{ln}({\bm k})$ is given by
([\[eq:berry-connection-matrix\]](#eq:berry-connection-matrix){reference-type="ref"
reference="eq:berry-connection-matrix"}) Finally, the matrix
$B_{nl}^{ac}$ has both orbital and spin contributions given by
$$\label{eq:B-ac-orb}
B_{nl}^{ac\,({\rm orb})}=
  \langle u_n\vert(\partial_aH)\vert\partial_c u_l\rangle
 -\langle\partial_c u_n\vert(\partial_aH)\vert u_l\rangle$$ and
$$\label{eq:B-ac-spin}
B_{nl}^{ac\,({\rm spin})}=-\frac{i\hbar^2}{m_e}\epsilon_{abc}
\langle u_n\vert\sigma_b\vert u_l\rangle.$$ The spin matrix elements
contribute less than 0.5% of the total $\rho_0^{\rm inter}$ of Te.
Expanding $H=\sum_m \vert u_m\rangle E_m \langle u_m\vert$ we obtain for
the orbital matrix elements
$$B_{nl}^{ac\,({\rm orb})}=-i\partial_a(E_n+E_l)A_{nl}^c \sum_m \Bigl\{ (E_n-E_m) A_{nm}^aA_{ml}^c -(E_l-E_m) A_{nm}^cA_{ml}^a \Bigr\}.
\label{eq:Bnl-sum}$$ This reduces the calculation of $B^{\text{(orb)}}$
to the evaluation of band gradients and off-diagonal elements of the
Berry connection matrix. Both operations can be carried out efficiently
in a Wannier-function basis following Ref. [@yates-prb07].

### `gyrotropic_task=-spin`: compute also the spin component of NOA and KME 

Unless this task is specified, only the orbital contributions are
calcuated in NOA and KME, thus contributions from
Eqs. ([\[eq:m-spin\]](#eq:m-spin){reference-type="ref"
reference="eq:m-spin"})
and ([\[eq:B-ac-spin\]](#eq:B-ac-spin){reference-type="ref"
reference="eq:B-ac-spin"}) are omitted.

## Electronic transport calculations with the `BoltzWann` module {#ch:boltzwann}

By setting $\verb#boltzwann#=\verb#TRUE#$, `postw90` will call the
`BoltzWann` routines to calculate some transport coefficients using the
Boltzmann transport equation in the relaxation time approximation.

In particular, the transport coefficients that are calculated are: the
electrical conductivity $\bm{\mathrm{\sigma}}$, the Seebeck coefficient
$\bm{\mathrm{S}}$ and the coefficient $\bm{\mathrm{K}}$ (defined below;
it is the main ingredient of the thermal conductivity).

The list of parameters of the `BoltzWann` module are summarized in
Table [\[parameter_keywords_bw\]](#parameter_keywords_bw){reference-type="ref"
reference="parameter_keywords_bw"}. An example of a Boltzmann transport
calculation can be found in the `wannier90` Tutorial.

**Note**: By default, the code assumes to be working with a 3D bulk
material, with periodicity along all three spatial directions. If you
are interested in studying 2D systems, set the correct value for the
`boltz_2d_dir` variable (see
Sec. [\[sec:boltz2ddir\]](#sec:boltz2ddir){reference-type="ref"
reference="sec:boltz2ddir"} for the documentation). This is important
for the evaluation of the Seebeck coefficient.

Please cite the following paper [@pizzi-cpc14] when publishing results
obtained using the `BoltzWann` module:

> G. Pizzi, D. Volja, B. Kozinsky, M. Fornari, and N. Marzari,\
> *BoltzWann: A code for the evaluation of thermoelectric and electronic
> transport properties with a maximally-localized Wannier functions
> basis*,\
> Comp. Phys. Comm. 185, 422 (2014), DOI:10.1016/j.cpc.2013.09.015.

### Theory {#sec:boltzwann-theory}

The theory of the electronic transport using the Boltzmann transport
equations can be found for instance in
Refs. [@ziman-book72; @grosso-book00; @mahan-itc06]. Here we briefly
summarize only the main results.

The current density $\bm{\mathrm{J}}$ and the heat current (or energy
flux density) $\bm{\mathrm{J}}_Q$ can be written, respectively, as
$$\begin{aligned}
  \bm{\mathrm{J}}   &= \bm{\mathrm{\sigma}}(\bm{\mathrm{E}} - \bm{\mathrm{S}} \bm{\mathrm{\nabla }}T) \\
  \bm{\mathrm{J}}_Q &= T \bm{\mathrm{\sigma }}\bm{\mathrm{S}} \bm{\mathrm{E}} - \bm{\mathrm{K}} \bm{\mathrm{\nabla }}T,
\end{aligned}$$ where the electrical conductivity
$\bm{\mathrm{\sigma}}$, the Seebeck coefficient $\bm{\mathrm{S}}$ and
$\bm{\mathrm{K}}$ are $3\times 3$ tensors, in general.

Note: the thermal conductivity $\bm{\mathrm{\kappa}}$ (actually, the
electronic part of the thermal conductivity), which is defined as the
heat current per unit of temperature gradient in open-circuit
experiments (i.e., with $\bm{\mathrm{J}}=0$) is not precisely
$\bm{\mathrm{K}}$, but
$\bm{\mathrm{\kappa }}= \bm{\mathrm{K}}-\bm{\mathrm{S}} \bm{\mathrm{\sigma }}\bm{\mathrm{S}} T$
(see for instance Eq. (7.89) of Ref. [@ziman-book72] or Eq. (XI-57b) of
Ref. [@grosso-book00]). The thermal conductivity $\bm{\mathrm{\kappa}}$
can be then calculated from the $\bm{\mathrm{\sigma}}$,
$\bm{\mathrm{S}}$ and $\bm{\mathrm{K}}$ tensors output by the code.

These quantities depend on the value of the chemical potential $\mu$ and
on the temperature $T$, and can be calculated as follows:
$$\begin{aligned}
_{ij}(\mu,T)&=e^2 \int_{-\infty}^{+\infty} d\varepsilon \left(-\frac {\partial f(\varepsilon,\mu,T)}{\partial \varepsilon}\right)\Sigma_{ij}(\varepsilon), \\
  [\bm{\mathrm{\sigma }}\bm{\mathrm{S}}]_{ij}(\mu,T)&=\frac e T \int_{-\infty}^{+\infty} d\varepsilon \left(-\frac {\partial f(\varepsilon,\mu,T)}{\partial \varepsilon}\right)(\varepsilon-\mu)\Sigma_{ij}(\varepsilon), \label{eq:boltz-sigmas}\\
  [\bm{\mathrm{K}}]_{ij}(\mu,T)&=\frac 1 T \int_{-\infty}^{+\infty} d\varepsilon \left(-\frac {\partial f(\varepsilon,\mu,T)}{\partial \varepsilon}\right)(\varepsilon-\mu)^2 \Sigma_{ij}(\varepsilon),\label{eq:boltz-thermcond}
\end{aligned}$$ where $[\bm{\mathrm{\sigma }}\bm{\mathrm{S}}]$ denotes
the product of the two tensors $\bm{\mathrm{\sigma}}$ and
$\bm{\mathrm{S}}$, $f(\varepsilon,\mu,T)$ is the usual Fermi--Dirac
distribution function
$$f(\varepsilon,\mu,T) = \frac{1}{e^{(\varepsilon-\mu)/K_B T}+1}$$ and
$\Sigma_{ij}(\varepsilon)$ is the Transport Distribution Function (TDF)
tensor, defined as
$$\Sigma_{ij}(\varepsilon) = \frac 1 V \sum_{n,\bm{\mathrm{k}}} v_i(n,\bm{\mathrm{k}}) v_j(n,\bm{\mathrm{k}}) \tau(n,\bm{\mathrm{k}}) \delta(\varepsilon - E_{n,k}).$$

In the above formula, the sum is over all bands $n$ and all states
$\bm{\mathrm{k}}$ (including spin, even if the spin index is not
explicitly written here). $E_{n,\bm{\mathrm{k}}}$ is the energy of the
$n-$th band at $\bm{\mathrm{k}}$, $v_i(n,\bm{\mathrm{k}})$ is the $i-$th
component of the band velocity at $(n,\bm{\mathrm{k}})$, $\delta$ is the
Dirac's delta function, $V = N_k \Omega_c$ is the total volume of the
system ($N_k$ and $\Omega_c$ being the number of $k$-points used to
sample the Brillouin zone and the unit cell volume, respectively), and
finally $\tau$ is the relaxation time. In the *relaxation-time
approximation* adopted here, $\tau$ is assumed as a constant, i.e., it
is independent of $n$ and $\bm{\mathrm{k}}$ and its value (in fs) is
read from the input variable `boltz_relax_time`.

### Files

#### `seedname_boltzdos.dat`

OUTPUT. Written by `postw90` if `boltz_calc_also_dos` is `true`. Note
that even if there are other general routines in `postw90` which
specifically calculate the DOS, it may be convenient to use the routines
in `BoltzWann` setting `boltz_calc_also_dos = true` if one must also
calculate the transport coefficients. In this way, the (time-demanding)
band interpolation on the $k$ mesh is performed only once, resulting in
a much shorter execution time.

The first lines are comments (starting with \# characters) which
describe the content of the file. Then, there is a line for each energy
$\varepsilon$ on the grid, containing a number of columns. The first
column is the energy $\varepsilon$. The following is the DOS at the
given energy $\varepsilon$. The DOS can either be calculated using the
adaptive smearing scheme[^6] if `boltz_dos_adpt_smr` is `true`, or using
a "standard" fixed smearing, whose type and value are defined by
`boltz_dos_smr_type` and `boltz_dos_smr_fixed_en_width`, respectively.
If spin decomposition is required (input flag `spin_decomp`), further
columns are printed, with the spin-up projection of the DOS, followed by
spin-down projection.

#### `seedname_tdf.dat`

OUTPUT. This file contains the Transport Distribution Function (TDF)
tensor $\bm{\mathrm{\Sigma}}$ on a grid of energies.

The first lines are comments (starting with \# characters) which
describe the content of the file. Then, there is a line for each energy
$\varepsilon$ on the grid, containing a number of columns. The first is
the energy $\varepsilon$, the followings are the components if
$\bm{\mathrm{\Sigma}}(\varepsilon)$ in the following order:
$\Sigma_{xx}$, $\Sigma_{xy}$, $\Sigma_{yy}$, $\Sigma_{xz}$,
$\Sigma_{yz}$, $\Sigma_{zz}$. If spin decomposition is required (input
flag `spin_decomp`), 12 further columns are provided, with the 6
components of $\bm{\mathrm{\Sigma}}$ for the spin up, followed by those
for the spin down.

The energy $\varepsilon$ is in eV, while $\bm{\mathrm{\Sigma}}$ is in
$\displaystyle\frac{1}{\hbar^2}\cdot\frac{\mathrm{eV}\cdot\mathrm{fs}}{\text{\AA}}$.

#### `seedname_elcond.dat`

OUTPUT. This file contains the electrical conductivity tensor
$\bm{\mathrm{\sigma}}$ on the grid of $T$ and $\mu$ points.

The first lines are comments (starting with \# characters) which
describe the content of the file. Then, there is a line for each
$(\mu,T)$ pair, containing 8 columns, which are respectively: $\mu$,
$T$, $\sigma_{xx}$, $\sigma_{xy}$, $\sigma_{yy}$, $\sigma_{xz}$,
$\sigma_{yz}$, $\sigma_{zz}$. (The tensor is symmetric).

The chemical potential is in eV, the temperature is in K, and the
components of the electrical conductivity tensor ar in SI units, i.e. in
1/$\Omega$/m.

#### `seedname_sigmas.dat`

OUTPUT. This file contains the tensor
$\bm{\mathrm{\sigma}}\bm{\mathrm{S}}$, i.e. the product of the
electrical conductivity tensor and of the Seebeck coefficient as defined
by Eq. [\[eq:boltz-sigmas\]](#eq:boltz-sigmas){reference-type="eqref"
reference="eq:boltz-sigmas"}, on the grid of $T$ and $\mu$ points.

The first lines are comments (starting with \# characters) which
describe the content of the file. Then, there is a line for each
$(\mu,T)$ pair, containing 8 columns, which are respectively: $\mu$,
$T$, $(\sigma S)_{xx}$, $(\sigma S)_{xy}$, $(\sigma S)_{yy}$,
$(\sigma S)_{xz}$, $(\sigma S)_{yz}$, $(\sigma S)_{zz}$. (The tensor is
symmetric).

The chemical potential is in eV, the temperature is in K, and the
components of the tensor ar in SI units, i.e. in A/m/K.

#### `seedname_seebeck.dat`

OUTPUT. This file contains the Seebeck tensor $\bm{\mathrm{S}}$ on the
grid of $T$ and $\mu$ points.

Note that in the code the Seebeck coefficient is defined as zero when
the determinant of the electrical conductivity $\bm{\mathrm{\sigma}}$ is
zero. If there is at least one $(\mu, T)$ pair for which
$\det \bm{\mathrm{\sigma}}=0$, a warning is issued on the output file.

The first lines are comments (starting with \# characters) which
describe the content of the file. Then, there is a line for each
$(\mu,T)$ pair, containing 11 columns, which are respectively: $\mu$,
$T$, $S_{xx}$, $S_{xy}$, $S_{xz}$, $S_{yx}$, $S_{yy}$, $S_{yz}$,
$S_{zx}$, $S_{zy}$, $S_{zz}$.

NOTE: therefore, the format of the columns of this file is different
from the other three files (elcond, sigmas and kappa)!

The chemical potential is in eV, the temperature is in K, and the
components of the Seebeck tensor ar in SI units, i.e. in V/K.

#### `seedname_kappa.dat`

OUTPUT. This file contains the tensor $\bm{\mathrm{K}}$ defined in
Sec. [3.1](#sec:boltzwann-theory){reference-type="ref"
reference="sec:boltzwann-theory"} on the grid of $T$ and $\mu$ points.

The first lines are comments (starting with \# characters) which
describe the content of the file. Then, there is a line for each
$(\mu,T)$ pair, containing 8 columns, which are respectively: $\mu$,
$T$, $K_{xx}$, $K_{xy}$, $K_{yy}$, $K_{xz}$, $K_{yz}$, $K_{zz}$. (The
tensor is symmetric).

The chemical potential is in eV, the temperature is in K, and the
components of the $\bm{\mathrm{K}}$ tensor are the SI units for the
thermal conductivity, i.e. in W/m/K.

## Generic Band interpolation {#ch:geninterp}

By setting $\verb#geninterp#=\verb#TRUE#$, `postw90` will calculate the
band energies (and possibly the band derivatives, if also
`geninterp_alsofirstder` is set to `TRUE`) on a generic list of $k$
points provided by the user.

The list of parameters of the Generic Band Interpolation module are
summarized in
Table [\[parameter_keywords_geninterp\]](#parameter_keywords_geninterp){reference-type="ref"
reference="parameter_keywords_geninterp"}. The list of input $k$ points
for which the band have to be calculated is read from the file named
`seedname_geninterp.kpt`. The format of this file is described below.

### Files

#### `seedname_geninterp.kpt`

INPUT. Read by `postw90` if `geninterp` is `true`.

The first line is a comment (its maximum allowed length is 500
characters).

The second line must contain `crystal` (or `frac`) if the $k$-point
coordinates are given in crystallographic units, i.e., in fractional
units with respect to the primitive reciprocal lattice vectors.
Otherwise, it must contain `cart` (or `abs`) if instead the $k-$point
coordinates are given in absolute coordinates (in units of 1/Å) along
the $k_x$, $k_y$ and $k_z$ axes.

*Note on units*: In the case of absolute coordinates, if $a_{lat}$ is
the lattice constant expressed in angstrom, and you want to represent
for instance the point $X=\frac {2\pi}{a_{lat}} [0.5, 0, 0]$, then you
have to input for its $x$ coordinate $k_x = 0.5 * 2 * \pi / a_{lat}$. As
a practical example, if $a_{lat}=4$Å, then $k_x = 0.78539816339745$ in
absolute coordinates in units of 1/Å.

The third line must contain the number $n$ of following $k$ points.

The following $n$ lines must contain the list of $k$ points in the
format

    kpointidx k1 k2 k3

where `kpointidx` is an integer identifying the given $k$ point, and
`k1`, `k2` and `k3` are the three coordinates of the $k$ points in the
chosen units.

#### `seedname_geninterp.dat` or ` seedname_geninterp_NNNNN.dat` {#sec:seedname.geninterp.dat}

OUTPUT. This file/these files contain the interpolated band energies
(and also the band velocities if the input flag `geninterp_alsofirstder`
is `true`).

If the flag `geninterp_single_file` is `true`, then a single file
`seedname_geninterp.dat` is written by the code at the end of the
calculation. If instead one sets `geninterp_single_file` to `false`,
each process writes its own output file, named
`seedname_geninterp_00000.dat`, ` seedname_geninterp_00001.dat`, ...

This flag is useful when one wants to parallelize the calculation on
many nodes, and it should be used especially for systems with a small
number of Wannier functions, when one wants to compute the bands on a
large number of $k$ points (if the flag `geninterp_single_file` is
`true`, instead, all the I/O is made by the root node, which is a
significant bottleneck).

**Important!** The files are not deleted before the start of a
calculation, but only the relevant files are overwritten. Therefore, if
one first performs a calculation and then a second one with a smaller
number of processors, care is needed to avoid to mix the results of the
older calculations with those of the new one. In case of doubt, either
check the date stamp in the first line of the
` seedname_geninterp_*.dat` files, or simply delete the
` seedname_geninterp_*.dat` files before starting the new calculation.

To join the files, on can simply use the following command:

    cat seedname_geninterp_*.dat > seedname_geninterp.dat

or, if one wants to remove the comment lines:

    rm seedname_geninterp.dat
    for i in seedname_geninterp_*.dat ; do grep -v \# "$i" >> \
    seedname_geninterp.dat ; done

The first few lines of each files are comments (starting with #),
containing a datestamp, the comment line as it is read from the input
file, and a header. The following lines contain the band energies (and
derivatives) for each band and $k$ point (the energy index runs faster
than the $k$-point index). For each of these lines, the first four
columns contain the $k$-point index as provided in the input, and the
$k$ coordinates (always in absolute coordinates, in units of 1/Å). The
fifth column contains the band energy.

If `geninterp_alsofirstder` is `true`, three further columns are
printed, containing the three first derivatives of the bands along the
$k_x$, $k_y$ and $k_z$ directions (in units of eV$\cdot$Å).

The $k$ point coordinates are in units of 1/Å, the band energy is in eV.

# Appendices

## Utilities {#ch:utilities}

The `wannier90` code is shipped with a few utility programs that may be
useful in some occasions. In this chapter, we describe their use.

### $\tt{kmesh.pl}$ {#sec:kmesh}

The `wannier90` code requires the definition of a full Monkhorst--Pack
grid of $k$ points. In the input file the size of this mesh is given by
means of the `mp_grid` variable. E.g., setting

    mp_grid = 4 4 4

tells `wannier90` that we want to use a $4\times 4\times 4$ $k$ grid.

One has then to specify (inside the `kpoints` block in the the
`seedname.win` file) the list of $k$ points of the grid. Here, the
`kmesh.pl` Perl script becomes useful, being able to generate the
required list.

The script can be be found in the `utility` directory of the
`wannier90` distribution. To use it, simply type:

    ./kmesh.pl nx ny nz

where `nx`, `ny` and `nz` define the size of the Monkhorst--Pack grid
that we want to use (for instance, in the above example of the
$4\times 4\times 4$ $k$ grid, `nx`$=$`ny`$=$`nz`$=$`<!-- -->`{=html}4).

This produces on output the list of $k$ points in Quantum Espresso
format, where (apart from a header) the first three columns of each line
are the $k$ coordinates, and the fourth column is the weight of each $k$
point. This list can be used to create the input file for the ab-initio
`nscf` calculation.

If one wants instead to generate the list of the $k$ coordinates without
the weight (in order to copy and paste the output inside the
`seedname.win` file), one simply has to provide a fourth argument on the
command line. For instance, for a $4\times 4\times 4$ $k$ grid, use

    ./kmesh.pl 4 4 4 wannier

and then copy the output inside the in the `kpoints` block in the
`seedname.win` file.

We suggest to always use this utility to generate the $k$ grids. This
allows to provide the $k$ point coordinates with the accuracy required
by `wannier90`, and moreover it makes sure that the $k$ grid used in the
ab-initio code and in `wannier90` are the same.

### $\tt{w90chk2chk.x}$[]{#sec:w90chk2chk label="sec:w90chk2chk"}

During the calculation of the Wannier functions, `wannier90` produces a
`.chk` file that contains some information to restart the calculation.

This file is also required by the `postw90` code. In particular, the
`postw90` code requires at least the `.chk` file, the `.win` input file,
and (almost always) the `.eig` file. Specific modules may require
further files: see the documentation of each module.

However, the `.chk` file is written in a machine-dependent format. If
one wants to run `wannier90` on a machine, and then continue the
calculation with `postw90` on a different machine (or with
`postw90` compiled with a different compiler), the file has to be
converted first in a machine-independent "formatted" format on the first
machine, and then converted back on the second machine.

To this aim, use the `w90chk2chk.x` executable. Note that this
executable is not compiled by default: you can obtain it by executing

    make w90chk2chk

in the main `wannier90` directory.

A typical use is the following:

1.  Calculate the Wannier functions with `wannier90`

2.  At the end of the calculation you will find a `seedname.chk` file.
    Run (in the folder with this file) the command

        w90chk2chk.x -export seedname

    or equivalently

        w90chk2chk.x -u2f seedname

    (replacing `seedname` with the seedname of your calculation).

    This command reads the `seedname.chk` file and creates a formatted
    file `seedname.chk.fmt` that is safe to be transferred between
    different machines.

3.  Copy the `seedname.chk.fmt` file (together with the `seedname.win`
    and `seedname.eig` files) on the machine on which you want to run
    `postw90`.

4.  On this second machine (after having compiled `w90chk2chk.x`) run

        w90chk2chk.x -import seedname

    or equivalently

        w90chk2chk.x -f2u seedname

    This command reads the `seedname.chk.fmt` file and creates an
    unformatted file `seedname.chk` ready to be used by `postw90`.

5.  Run the `postw90` code.

### $\tt{PL\_assessment}$ {#sec:pl_assessment}

The function of this utility is to assess the length of a principal
layer (in the context of a Landauer-Buttiker quantum conductance
calculation) of a periodic system using a calculation on a single unit
cell with a dense k-point mesh.

The utility requires the real-space Hamiltonian in the MLWF basis,
`seedname_hr.dat`.

The `seedname_hr.dat` file should be copied to a directory containing
executable for the utility. Within that directory, run:

    \$> ./PL_assess.x  nk1 nk2 nk3 num_wann 

where:

`nk1` is the number of k-points in x-direction `nk2` is the number of
k-points in y-direction `nk3` is the number of k-points in z-direction
`num_wann` is the number of wannier functions of your system

e.g.,

    \$> ./PL_assess.x  1 1 20 16

Note that the current implementation only allows for a single k-point in
the direction transverse to the transport direction.

When prompted, enter the seedname.

The programme will return an output file `seedname_pl.dat`, containing
four columns

1.  Unit cell number, $R$

2.  Average 'on-site' matrix element between MLWFs in the home unit
    cell, and the unit cell $R$ lattice vectors away

3.  Standard devaition of the quantity in (2)

4.  Maximum absolute value in (2)

### $\tt{w90vdw}$ {#sec:w90vdw}

This utility provides an implementation of a method for calculating van
der Waals energies based on the idea of density decomposition via MLWFs.

For theoretical details, please see the following publication and
references therein:

Lampros Andrinopoulos, Nicholas D. M. Hine and Arash A. Mostofi,
"Calculating dispersion interactions using maximally localized Wannier
functions", *J. Chem. Phys.* **135**, 154105 (2011).

For further details of this program, please see the documentation in
`utility/w90vdw/doc/` and the related examples in
`utility/w90vdw/examples/`.

### $\tt{w90pov}$ {#sec:w90pov}

An utility to create Pov files (to render the Wannier functions using
the PovRay utility) is provided inside `utility/w90pov`.

Please refer to the documentation inside `utility/w90pov/doc` for more
information.

### $\tt{k\_mapper.py}$ {#sec:k_mapper}

The `wannier90` code requires the definition of a full Monkhorst--Pack
grid of $\mathbf{k}$-vectors, which can be obtained by means of the
`kmesh.pl` utility. In order to perform a GW calculation with the Yambo
code, you need to perform a nscf calculation on a grid in the
irreducible BZ. Moreover, you may need a finer grid to converge the GW
calculation than what you need to interpolate the band structure. The
`k_mapper.py` tools helps in finding the $\mathbf{k}$-vectors indexes of
a full grid needed for interpolation into the reduced grid needed for
the GW calculation with Yambo. Usage:

`path/k_mapper.py nx ny nz QE_nscf_output`

where `path` is the path of `utility` folder, `QE_nscf_output` is the
path of the QE nscf output file given to Yambo.

### $\tt{gw2wannier90.py}$

This utility allows to sort in energy the input data of `wannier90`
(e.g. overlap matrices and energy eigenvalues). `gw2wannier90.py` allows
to use `wannier90` at the $G_0W_0$ level, where quasi-particle
corrections can change the energy ordering of eigenvalues (Some
`wannier90` modules require states to be ordered in energy). Usage:

`path/gw2wannier90.py seedname options` where `path` is the path of
`utility` folder.

Available options are:

    mmn, amn, spn, unk, uhu, uiu,
    spn_formatted, unk_formatted, uhu_formatted, uiu_formatted,
    write_formatted

If no options are specified, all the files
(`mmn, amn, spn, UNK, uHu, uIu`) are considered.

Binary (unformatted Fortran) files are supported, though not
reccommended, since they are compiler-dependent. A safer choice is to
use (bigger) formatted files, with options:

`spn_formatted, uiu_formatted, uhu_formatted, unk_formatted`

In default, the output format is the same as the input format. To
generate formatted files with unformatted input, use option:
`write_formatted` []{#sec:w90aaa label="sec:w90aaa"}

### $\tt{w90spn2spn.x}$[]{#sec:w90spn2spn label="sec:w90spn2spn"}

The interface between ab-initio code and `wannier90` (e.g.
`pw2wannier90.x`) can produce a `.spn` file that is used by `postw90` to
calculate some spin related quantities.

The `.spn` file can be written in a machine-dependent or a
machine-independent format depending on the input parameter
`spn_formatted` (the default is `false` which means the `.spn` file is
machine-dependent) of the `pw2wannier90.x`. If a `.spn` file has been
generated on a machine with machine-dependent format, and then one wants
to continue the calculation with `postw90` on a different machine (or
with `postw90` compiled with a different compiler), the file has to be
converted first in a machine-independent "formatted" format on the first
machine.

To this aim, use the `w90spn2spn.x` executable. Note that this
executable is not compiled by default: you can obtain it by executing

    make w90spn2spn

in the main `wannier90` directory.

A typical use is the following:

1.  Calculate the `.spn` file, e.g. by `pw2wannier90.x`

2.  At the end of the calculation you will find a `seedname.spn` file.
    If the file is "unformatted", run (in the folder with this file) the
    command

        	w90spn2spn.x -export seedname
        	

    or equivalently

        	w90spn2spn.x -u2f seedname
        	

    (replacing `seedname` with the seedname of your calculation).

    This command reads the `seedname.spn` file and creates a formatted
    file `seedname.spn.fmt` that is safe to be transferred between
    different machines.

3.  Copy the `seedname.spn.fmt` file on the machine on which you want to
    run `postw90`.

4.  On this second machine (after having compiled `w90spn2spn.x`) run

        	w90spn2spn.x -import seedname
        	

    or equivalently

        	w90spn2spn.x -f2u seedname
        	

    This command reads the `seedname.spn.fmt` file and creates an
    unformatted file `seedname.spn` ready to be used by `postw90`.

5.  Run the `postw90` code.

Note if `spn_formatted` is set to `true` in both `pw2wannier90.x` and
`postw90` input files, then the `.spn` file will be written and read as
"formatted", so `w90spn2spn.x` is not needed. However, if an
"unformatted" `seedname.spn` has been created and you do not want to
rerun `pw2wannier90.x`, then `w90spn2spn.x` can be useful. Also, once a
"formatted" `seedname.spn` has been generated, the step 4 can be skipped
if `spn_formatted` is set to `true` in `postw90` input file
`seedname.win`.

### `write_pdwf_projectors.py`

A python script to extract projectors from a `UPF` file and write them
into a `pw2wannier90.x`-recognizable `dat` file, which can be used to
compute `amn` using pseudo-atomic orbital projection.

Usage:

`path/write_pdwf_projectors.py UPF_filename`

where `path` is the path of `utility` folder, `UPF_filename` is the path
of a `UPF` file.

The script serves as a reference for writing the `dat` file, you can
generate your own pseudo-atomic orbitals by some other codes and use
them to compute `amn`.

## Frequently Asked Questions {#chap:faq}

### General Questions

#### What is `wannier90`?

`wannier90` is a computer package, written in Fortran90, for obtaining
maximally-localised Wannier functions, using them to calculate
bandstructures, Fermi surfaces, dielectric properties, sparse
Hamiltonians and many things besides.

#### Where can I get `wannier90`?

The most recent release of `wannier90` is always available on our
website <http://www.wannier.org>.

#### Where can I get the most recent information about `wannier90`?

The latest news about `wannier90` can be followed on our website
<http://www.wannier.org>.

#### Is `wannier90` free?

Yes! `wannier90` is available for use free-of-charge under the GNU
General Public Licence. See the file `LICENSE` that comes with the
`wannier90` distribution or the GNU hopepage at <http://www.gnu.org>.

### Getting Help

#### Is there a Tutorial available for `wannier90`?

Yes! The `examples` directory of the `wannier90` distribution contains
input files for a number of tutorial calculations. The ` doc` directory
contains the accompanying tutorial handout.

#### Where do I get support for `wannier90`?

There are a number of options:

1.  The `wannier90` User Guide, available in the `doc` directory of the
    distribution, and from the webpage
    (<http://www.wannier.org/user_guide.html>)

2.  The `wannier90` webpage for the most recent announcements
    (<http://www.wannier.org>)

3.  The `wannier90` mailing list (see
    <http://www.wannier.org/forum.html>)

#### Is there a mailing list for `wannier90`?

Yes! You need to register: go to <http://www.wannier.org/forum.html> and
follow the instructions.

### Providing Help: Finding and Reporting Bugs

#### I think I found a bug. How do I report it?

-   Check and double-check. Make sure it's a bug.

-   Check that it is a bug in `wannier90` and not a bug in the software
    interfaced to `wannier90`.

-   Check that you're using the latest version of `wannier90`.

-   Send us an email. Make sure to describe the problem and to attach
    all input and output files relating to the problem that you have
    found.

#### I have got an idea! How do I report a wish?

We're always happy to listen to suggestions. Email your idea to the
`wannier90` developers.

#### I want to help! How can I contribute to `wannier90`?

Great! There's always plenty of functionality to add. Email us to let us
know about the functionality you'd like to contribute.

#### I like `wannier90`! Should I donate anything to its authors?

Our Swiss bank account number is\... just kidding! There is no need to
donate anything, please just cite our paper in any publications that
arise from your use of `wannier90`:

> \[ref\] G. Pizzi, V. Vitale, R. Arita, S. Blügel, F. Freimuth, G.
> Géranton, M. Gibertini, D. Gresch, C. Johnson, T. Koretsune, J.
> Ibañez-Azpiroz, H. Lee, J.M. Lihm, D. Marchand, A. Marrazzo, Y.
> Mokrousov, J.I. Mustafa, Y. Nohara, Y. Nomura, L. Paulatto, S. Poncé,
> T. Ponweiser, J. Qiao, F. Thöle, S.S. Tsirkin, M. Wierzbowska, N.
> Marzari, D. Vanderbilt, I. Souza, A.A. Mostofi, J.R. Yates,\
> Wannier90 as a community code: new features and applications, *J.
> Phys. Cond. Matt.* **32**, 165902 (2020)\
> <https://doi.org/10.1088/1361-648X/ab51ff>

If you are using versions 2.x of the code, cite instead:

> \[ref\] A. A. Mostofi, J. R. Yates, G. Pizzi, Y.-S. Lee, I. Souza,
> D. Vanderbilt and N. Marzari,\
> An updated version of `wannier90`: A Tool for Obtaining
> Maximally-Localised Wannier Functions, *Comput. Phys. Commun.*
> **185**, 2309 (2014)\
> <http://doi.org/10.1016/j.cpc.2014.05.003>

### Installation

#### How do I install `wannier90`?[]{#sec:installation label="sec:installation"}

Follow the instructions in the file `README.install` in the main
directory of the `wannier90` distribution.

#### Are there `wannier90` binaries available?

Not at present.

#### Is there anything else I need?

Yes. `wannier90` works on top of an electronic structure calculation.

At the time of writing there are public, fully functioning, interfaces
between `wannier90` and [pwscf]{.smallcaps}, abinit
(<http://www.abinit.org>), siesta (<http://www.icmab.es/siesta/>), VASP
(<https://www.vasp.at>), Wien2k (<http://www.wien2k.at>), fleur
(<http://www.fleur.de>), OpenMX (<http://www.openmx-square.org/>), GPAW
(<https://wiki.fysik.dtu.dk/gpaw/>).

To use `wannier90` in combination with [pwscf]{.smallcaps} code (a
plane-wave, pseudopotential, density-functional theory code, which is
part of the `quantum-espresso` package) you will need to download
[pwscf]{.smallcaps} from the webpage <http://www.quantum-espresso.org>.
Then compile [pwscf]{.smallcaps} and the `wannier90` interface program
`pw2wannier90`. For instructions, please refer to the documentation that
comes with the `quantum-espresso` distribution.

For examples of how to use [pwscf]{.smallcaps} and `wannier90` in
conjunction with each other, see the `wannier90` Tutorial.

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
