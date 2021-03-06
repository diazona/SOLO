SOLO: Saturation physics at One Loop Order
==================================================

This is the program used to calculate the complete next-to-leading cross section
for inclusive hadron production in pA collisions, described in the paper

Anna M. Stasto, Bo-Wen Xiao, David Zaslavsky  
"Towards the Test of Saturation Physics Beyond Leading Logarithm"  
Phys. Rev. Lett. 112, 012302 (2014)  
arXiv:1307.4057 [hep-ph]

Please cite this paper if you use the results of the code in a publication.

Official documentation for the program is kept at
https://diazona.github.io/SOLO/.

Maintenance Status
==================================================

SOLO is written and maintained by me, David Zaslavsky. I no longer work
in academia, so the code is not actively being developed, but I'm happy
to answer questions about the program to the best of my ability. You
can email me with questions about SOLO at diazona@ellipsix.net or
the email address listed in my GitHub profile, if different.

Installation
==================================================

The quickest and intended way to compile the code is as follows: first, ensure
that [git](https://git-scm.com/), [GSL](https://www.gnu.org/software/gsl/),
[MuParser](https://github.com/beltoforion/muparser), and
[CMake](https://cmake.org/) are properly installed (as well as a C++ compiler
such as [GCC](https://www.gnu.org/software/gcc/)). Then, in a shell, run
the following commands:

    git clone https://github.com/diazona/SOLO.git
    cd SOLO
    git submodule init
    git submodule update

Then, download the MSTW PDF code from http://mstwpdf.hepforge.org/code/code.html.
Extract the files `mstwpdf.cc` and `mstwpdf.h` from the tarball and place
them in SOLO's `src/` directory. (We are not authorized to distribute
the MSTW PDF interface as part of SOLO, which is why this has to be done
manually.) Then, from the directory `SOLO/` (the parent
directory of `src/`), run these commands:
    
    mkdir build
    cd build
    cmake .. && make

At the end of this you should have a `build/src/` directory containing the
program `oneloopcalc` and other programs.

Other ways of installing the program (e.g. without git, or without Cmake)
are described in [the documentation](https://diazona.github.io/SOLO/).

Running the program
==================================================

## Additional files

In order to run the program, you will need two additional files:

- The grid file for the MSTW 2008 PDF at NLO, from the paper
    
    A. D. Martin, W. J. Stirling, R. S. Thorne and G. Watt,  
    "Parton distributions for the LHC",  
    Eur. Phys. J. C 63 (2009) 189-285  
    arXiv:0901.0002 [hep-ph].

The filename is `mstw2008nlo.00.dat`, and it can be downloaded as part of an
archive at the MSTW PDF site http://mstwpdf.hepforge.org/code/code.html.

- The data file for the DSS fragmentation functions at NLO, from the paper

    Daniel de Florian, Rodolfo Sassot, Marco Stratmann  
    "Global analysis of fragmentation functions for pions and kaons and their uncertainties"  
    Phys. Rev. D 75, 114010 (2007)  
    arXiv:hep-ph/0703242

The filename is `PINLO.DAT`. Unfortunately we (authors of SOLO) are not aware of
a website where this file is directly available.

## Invocation

The program is invoked as

    oneloopcalc <options>

where the <options> can include any number of the following, in any order:

- Hard factor group specifications

    These tell the program which terms to calculate. A hard factor _group_
    specification is made of any number of individual hard factor specifications
    separated by commas, for example

        m.h02qq,m.h02gg

    The program will calculate the results for all the terms in the group and
    display a total for each group. You can name a group by prefixing the
    specification with a colon, like

        lo:m.h02qq,m.h02gg

    The name will be used to label a column in the results table printed
    when the program finishes.
    
    An _individual_ hard factor specification is a string like "p.h02qq" or
    "m.h16gg". The "p." at the beginning specifies the position space
    implementation, "r." specifies a position space implementation with the
    angular integral already done, and the "m." specifies the momentum space
    implementation. The prefix can be omitted, in which case position space
    is taken as the default. (Not recommended, as position space is highly
    inaccurate for some terms.) For compatibility with older versions, the
    program also accepts a colon instead of the period (like "p:h02qq").
    
    The rest of the string gives the name of a hard factor. The canonical set
    of possible names that can be used with a "p." prefix are all the return
    values from the get_name() methods in hardfactors_position.h, and similarly
    for "r." with hardfactors_radial.h and "m." with hardfactors_momentum.h.
    Here's a near-complete list:
    
        p.h02qq m.h02qq p.h12qq r.h12qq p.h14qq m.h14qq
        p.h02gg m.h02gg p.h12gg r.h12gg p.h12qqbar m.h12qqbar p.h16gg m.h16gg
        p.h112qg r.h112qg p.h122qg r.h122qg p.h14qg m.h14qg
        p.h112gq r.h112gq p.h122gq r.h122gq p.h14gq m.h14gq
        
    and also these, which go beyond what is in the paper:
        
        r.h12qq.1 r.h12qq.1A r.h12qq.1B r.h12qq.2 r.h12qq.3 r.h012qqexp
        m.h1qqexact m.h1ggexact
        
    The names are case-insensitive.

    It's also possible to specify the group of all leading order terms using
    the shortcut "lo", which is equivalent to
    
        m.h02qq,m.h02gg
        
    or the group of all next-to-leading order terms using the shortcut
    "nlo.std", which is equivalent to
    
        r.h12qq,m.h14qq,r.h12gg,m.h12qqbar,m.h16gg,r.h112gq,r.h122gq,m.h14gq,
        r.h112qg,r.h122qg,m.h14qg
    
    or the group of all next-to-leading order terms using the high-pT expansion
    for the diagonal channels using the shortcut "nlo.hipt", which is
    equivalent to

        m.h1qqexact,m.h1ggexact,r.h112gq,r.h122gq,m.h14gq,r.h112qg,r.h122qg,
        m.h14qg

    The shortcut "nlo" will choose between these latter two options
    automatically: "nlo.std" if approximate kinematics are in use
    (exact_kinematics = 0 in the Context), or "nlo.hipt" if exact kinematics
    are in use (exact_kinematics = 1).

    These shortcuts are defined in oneloopcalc.cpp. The default if no hard
    factor groups are specified on the command line is
    
        lo nlo

- Configuration file names

    Configuration files contain parameters for the program, in the format
    
        key1=value1
        key2=value2
        
    and so on. Keys are case-insensitive. The canonical list of keys which
    are used is the code in `context.cpp`. Here's a mostly-complete list:
    
    - `A` (no default)

        the mass number

    - `abserr` (default 1e-20)

        the absolute error at which to stop an integration, for strategies which
        use this termination condition

    - `alphas` (default 0.2)

        value for the fixed coupling

    - `beta` (default `11 - 2*Nf/3`)

        coefficient for the LO running coupling

    - `c` (no default)

        the centrality coefficient, 0-1

    - `c0r_optimization` (default true)

        if the factorization scale scheme is c0r, whether to skip calculating
        terms that should be zero

    - `CF` (default 1.5)

        the color factor

    - `coupling_type` (default fixed)

        "fixed" or "running"

    - `cubature_iterations` (default 1000000)

        number of calls to use for cubature integration

    - `exact_kinematics` (default false)

        whether to use exact kinematic expressions

    - `factorization_scale` (default fixed)

        "`fixed`" or "`4pT2`" or "`CpT2`" or "`c0r`" to specify how to set the
        factorization scale

    - `factorization_scale_coefficient` (no default)

        if factorization_scale is "CpT2", this is the coefficient to multiply by
        pT2 to get mu2

    - `ff_filename` (default `PINLO.DAT`)

        filename to read DSS FF data from

    - `gammaMV` (default 1)

        the anomalous dimension in the MV gluon distribution

    - `gdist_momentum_filename` (no default)

        file to read the momentum data for a gluon distribution from

    - `gdist_position_filename` (no default)

        file to read the position data for a gluon distribution from

    - `gdist_subinterval_limit` (default 10000)

        number of subdivisions to use when integrating a position gluon
        distribution

    - `gdist_type` (default `GBW`)

        the type of the gluon distribution, "`GBW`", "`MV`", "`fMV`", "`file`", or
        "`gbw+file`"

    - `hadron` (no default)

        the type of hadron detected, "`pi-`", "`pi0`", or "`pi+`"

    - `inf` (default 40)

        the cutoff used for integration over a theoretically infinite region

    - `integration_strategy` (default `VEGAS`)

        the integration type to use, "`MISER`", "`VEGAS`" (best), or "`QUASI`"

    - `lambda` (default 0.288)

        the exponent in the definition of the saturation scale

    - `lambdaMV` (default 0.241)

        the parameter in the MV gluon distribution, in GeV

    - `lambdaQCD` (default 0.2428711 = sqrt(0.0588))

        QCD lambda in GeV, used in the running coupling

    - `miser_iterations` (default 10000000 = 1e7)

        number of iterations to use in MISER integration

    - `mu2` (default 10)

        factorization scale in GeV, if factorization_scale is "fixed"

    - `Nc` (default 3)

        number of colors

    - `Nf` (default 3)

        number of flavors

    - `pdf_filename` (default `mstw2008nlo.00.dat`)

        filename to read MSTW PDF from

    - `projectile` (no default)

        the type of projectile, "deuteron" or "proton"

    - `pseudorandom_generator_seed` (default 0)

        seed for the GSL random number generator

    - `pseudorandom_generator_type` (default `mt19937`)

        algorithm to use for generating random numbers; allowed values are in
        the GSL documentation

    - `pT` (no default)

        comma-separated list of transverse momenta

    - `quasirandom_generator_type` (default `halten`)

        algorithm to use for generating quasirandom numbers for QMC integration;
        allowed values are in the GSL documentation

    - `quasi_iterations` (default 1000000)

        the number of iterations at which to stop quasi Monte Carlo integration

    - `regulator` (default 1)

        the position of the Landau pole for the regulated LO running coupling

    - `relerr` (default 0)

        the relative error at which to stop an integration, for strategies which
        use this termination condition

    - `satscale_source` (default `extract from momentum`)

        for a file gluon distribution, how to extract the saturation scale;
        allowed values are "`analytic`" (Q0²(x0/x)^λ), "`extract from momentum`"
        which determines the saturation scale by finding the momentum where the
        gluon distribution equals a fixed fraction of its value at a reference
        momentum, and "`extract from position`" which finds the radius where the
        gluon distribution equals a fixed threshold value

    - `satscale_threshold` (no default)

        if satscale_source is "`extract from momentum`" or
        "`extract from position`", this is the fixed threshold value (or fraction
        of its value at a reference point, in the momentum case) that the gluon
        distribution should equal at the saturation scale

    - `Sperp` (default 1)

        cross-sectional area of the hadron

    - `sqs` (no default)

        sqrt(s), the collider's CM energy

    - `TR` (default 0.5)

        group coefficient

    - `vegas_incremental_iterations` (default 1000000)

        number of function evaluations to use in each step of the VEGAS Monte
        Carlo algorithm after the first

    - `vegas_initial_iterations` (default 100000)

        number of function evaluations to use to refine the grid in the first
        step of the VEGAS algorithm

    - `x0` (default 0.000304)

        the fit parameter from the definition of the saturation scale

    - `Y` (no default)

        comma-separated list of rapidities (in the center of mass frame) to run
        the calculation at

    The configuration files have to at least set `A`, `c`, `sqs`, and `Y`, and
    also `pT` if no transverse momenta are specified as command line arguments.
        
- Transverse momentum values

    Any numbers given as command line arguments are put together into one big
    list of transverse momentum values to run the calculation at. If a
    comma-separated list of numbers is given, then it will be split apart and
    each number added to the one big list. There's no significance to putting
    certain `pT` values together and others not. (`0.5 0.7 0.8,0.9` and
    `0.5,0.7 0.8 0.9` are exactly equivalent.) Any `pT` values specified on the
    command line will _replace_ `pT` values specified in the config file, if there
    is one.

- Literal options

    - `--separate`

        Print out the results for each individual hard factor, not just the total
        for each hard factor groups

    - `--minmax`

        Track and print out the minimum and maximum values of kinematic variables

    - `--trace-gdist`

		Print out parameters and values for every call to the gluon distribution.
		Output goes to the file trace_gdist.output in the working directory.
		(Expect this file to grow to several hundred megabytes.)

    - `--trace=var1,var2,...`

		Print out selected variables from the integration context after every
		single evaluation of the function. The output goes to the file `trace.output`
		in the working directory. (Expect this file to grow to several megabytes.)
		The allowable variables are those in `ictx_var_list.inc`, or you can use
		"`--trace=all`" or "`--trace=*`" to print out all available variables.

Structure
==================================================

Basically the code runs as follows:

1. Collect the command line options and settings from configuration files and
   put everything into a `ResultsCalculator`
2. For each combination of `pT` and `Y`, and for each hard factor group:

    1.  Create an `Integrator` with the current values of `pT` and `y` and the current
        hard factor group

    2. The Integrator calls the GSL Monte Carlo integration routine

    3. For each time the MC routine evaluates the function
    
        1. Update the variables in the `IntegrationContext`

        2. Go through the list of `HardFactor` instances in the current group
            and get a value from each one
        3. Return the total value

    d. Store the value and error bound returned from the Monte Carlo

3. Print out all the results


Files
==================================================

Source code for the program itself:

`oneloopcalc.cpp`  
Main program

`log.h`  
Declares an output stream to write status messages to

`gsl_exception.h`  
Declares an exception to be thrown when GSL reports an error

`hardfactors.h`  
A class that abstractly represents a hard factor (i.e. an expression to be
integrated)

`hardfactors_momentum.h`  
`hardfactors_momentum.cpp`  
Implementation of the momentum space hard factors (terms)

`hardfactors_position.h`  
`hardfactors_position.cpp`  
Implementation of the position space hard factors (terms)

`gluondist.h`  
`gluondist.cpp`  
Implementations of the gluon distributions

`gluondist_driver.cpp`  
A program to print out values from the gluon distributions

`coupling.h`  
`coupling.cpp`  
Implementations of the fixed and LO running couplings

`factorizationscale.h`  
`factorizationscale.cpp`  
Implementations of the various schemes for the factorization scale

`integrationcontext.h`  
`integrationcontext.cpp`  
A class that stores the kinematic variables used in the calculation. The
values stored in this get updated every time the function is evaluated.

`integrator.h`  
`integrator.cpp`  
A class that stores the parameters for the integral and actually calls the
GSL Monte Carlo integration functions

`integrationtype.h`  
`integrationtype.cpp`  
Definitions of integration types. An integration type specifies how many
dimensions are in the Monte Carlo integral and what the limits are.

`utils.h`  
`utils.cpp`  
Some string and list processing functions

`ictx_var_list.inc`  
Variables from the integration context, listed in a separate file as a
preprocessor hack of sorts

Source code and object code for other things used by the program:
    
`dss_pinlo.h`  
`dss_pinlo.cpp`  
A C++ interface to the DSS fragmentation functions

`dss_pinlo_test.cpp`  
Test program for the DSS FF interface

`interp2d.h`  
`libinterp2d.a`  
A 2D interpolation library compatible with the GSL
Full source code at https://github.com/diazona/interp2d

`quasimontecarlo.h`  
`libquasimontecarlo.a`  
A library for quasi Monte Carlo integration compatible with the GSL
Full source code at https://github.com/diazona/quasimontecarlo

Source code for other things used by the program, written by other people:
    
`mstwpdf.h`  
`mstwpdf.cc`  
A C++ interface to the MSTW PDFs

`cubature.h`  
`cubature.c`  
A library for multidimensional cubature (deterministic integration)
Not currently used

Non-source code files:
    
`CMakeFiles.txt`  
Instructions for the build system, CMake

`PINLO.DAT`  
DSS fragmentation function data

`pinlo_extended.dat`  
DSS fragmentation function data with extrapolation to lower z

`mstw2008nlo.00.dat`  
MSTW PDF data
