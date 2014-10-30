------------------------------------------------------------------------------------------------
Alpha expansion for Robust P^N potentials
------------------------------------------------------------------------------------------------

Authors: Pushmeet Kohli, Lubor Ladicky, Philip H.S.Torr
Oxford Brookes University


------------------------------------------------------------------------------------------------
Licence
------------------------------------------------------------------------------------------------

This software library implements the alpha expansion for robust P^N potentials described in

P. Kohli, L. Ladicky, and P. Torr. Graph cuts for minimizing robust higher order potentials.
Technical report, Oxford Brookes University, UK., 2008.

P. Kohli, L. Ladicky, and P. Torr. Robust higher order potentials for enforcing label
consistency. In CVPR, 2008.

------------------------------------------------------------------------------------------------

The library uses max-flow code described in

Yuri Boykov and Vladimir Kolmogorov. An Experimental Comparison of Min-Cut/Max-Flow Algorithms
for Energy Minimization in Vision. In IEEE Transactions on Pattern Analysis and Machine
Intelligence (PAMI), September 2004

------------------------------------------------------------------------------------------------

The alpha expansion algorithm is explained in

Y. Boykov, O. Veksler, and R. Zabih. Fast approximate energy minimization via graph cuts.
PAMI, 23(11):1222–1239, 2001.

------------------------------------------------------------------------------------------------

The code is free to use for research purposes. If you use this software for research purposes,
you should cite above papers in any resulting publication.

------------------------------------------------------------------------------------------------


------------------------------------------------------------------------------------------------
Files included
------------------------------------------------------------------------------------------------

README.txt        - this readme file

robustpn.cpp      - source file explaining how to use the library

energy.h          - header file with energy class definition

aexpand.h         - header file with alpha expansion code

graph.h, block.h  - header files containg max-flow code

------------------------------------------------------------------------------------------------


------------------------------------------------------------------------------------------------
Contact Information
------------------------------------------------------------------------------------------------

Email:
	lladicky@brookes.ac.uk  (Lubor Ladicky)
	pkohli@microsoft.com    (Pushmeet Kohli)

------------------------------------------------------------------------------------------------
