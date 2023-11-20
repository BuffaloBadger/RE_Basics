
Appendix I: when the corresponding SCoRE appendices are available, add links to the videos and to Matlab and Python functions for solving the example ATEs.



General

* add underlined variable symbols as vectors and bold symbols as matrices to nomenclature
* Revise preface once SCoRE is available.
* As SCoRE preparation assignment video is added for each example, add link to the numerical calculations videos (python) and (Matlab)
* Add link for reporting errors
* Revise any tables or text where scientific notation uses E to use x 10^x instead.
* Change code blocks that source fmt_tibble_col.R to use the copy in the Code folder.
* Double-check cross-references.
* In each chapter, describe the differences in the analytical procedures for the chapter problems under the examples header. (e.g.
in example 1 the reactor response is known and the input needed to cause that response is requested while in example 2, the input is 
known and the response is requested.)

Part 1 - Basic Information
* Change format of examples in chapters 1 through 3 to match chapter 5.

Part 2 - Rate Expressions
* Change format of examples in chapter 4 to match chapter 5.
* Revise example 4.4 to use linear least squares for finding Arrhenius parameters

Part 3 - Reactor Models

* Check the transient form of the energy balance on an exchange fluid that undergoes phase change.
* Get citation for transient momentum balance in Chapter 6
* Add an example of code for solving the reactor design equations to Chapter 6.

Part 4 - Kinetics Data Analysis
* Revise all chapters to refer to Chapter 7
* Create a standard formulation and calculations (see isolated reactor analysis) and re-write all examples using it.

Part 5 - Modeling Isolated Ideal Reactors

Part 6 - Modeling Ideal Reactor Systems

Part 7 - Modeling Non-Ideal Reactors

Appendices

* Check the transient form of the energy balance on an exchange fluid that undergoes phase change.
* Add illustration of laminar, turbulent and plug flow to prerequsite knowledge section on fluids; see slide 18 in VO Slides - RE Basics 1.3
* Fix figures in Appendix D so they all have the same size.
* Add linear least squares to prerequisite knowledge appendix.
* Add a section to Appendix G discussing packed beds and the psuedo-homogeneous assumption
* Change references to "problems" to references to "reaction engineering (or kinetics) assignments." 

Text cut from bottom of _quarto.yml

  pdf:
    documentclass: scrreprt
    include-in-header: 
      text: |
        \usepackage{cancel}
