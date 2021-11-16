#if 1
C     -- For completeness, we also include below the earlier routine
C     versions (with bug fixes).      
C
C     ##################################################################
C     Since circa 2010 SPARSPAK is no longer available through netlib.
C      
C     See original descriptions in George & Liu (1981),
C     Computer Solution of Large Sparse Positive-Definite Systems.
C     or (available as PDF) George, Lui & Ng (1994),
C     Computer Solution of Sparse Linear Systems.
C
C     ##################################################################
C
C     Here we supply routines from SPARSPAK that:
C     A. Carry out Reverse Cuthill--McKee node ordering (Ch 4 of book).
C        Main routine is RCM
C     B. Carry out nested dissection ordering           (Ch 8 of book).
C        Main routine is GENND
C
C     The subset below is taken from
C     https://people.sc.fsu.edu/~jburkardt/f77_src/sparsepak/sparsepak.f
C     which contains some later F77-style changes.  Quite a few errors
C     in FNDSEP (which were still in the 1994 book) were also fixed.
C
C     BUT! found an error in routine rcm - this is noted and fixed.
C
C     ##################################################################
      
      
      subroutine rootls ( root, xadj, adjncy, mask, nlvl, xls, ls )

c*********************************************************************72
c
cc ROOTLS generates the level structure rooted at a given node.
c
c  Discussion:
c
c    Only those nodes for which MASK is nonzero will be considered.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  root, the node at which the level structure is to be rooted.
c
c  (xadj, adjncy) - adjacency structure pair for the
c  given graph.
c
c  mask - is used to specify a section subgraph. nodes
c  with mask(i)=0 are ignored.
c
c  output parameters -
c  nlvl - is the number of levels in the level structure.
c  (xls, ls) - array pair for the rooted level structure.
c
      implicit none
      
      integer adjncy(*)
      integer ls(*), mask(*), xls(*)
      integer xadj(*), i, j, jstop, jstrt, lbegin
      integer ccsize, lvlend, lvsize, nbr, nlvl
      integer node, root
c
c  Initialization.
c
      mask(root) = 0
      ls(1) = root
      nlvl = 0
      lvlend = 0
      ccsize = 1
c
c  LBEGIN is the pointer to the beginning of the current
c  level, and lvlend points to the end of this level.
c
200   lbegin = lvlend + 1
      lvlend = ccsize
      nlvl = nlvl + 1
      xls(nlvl) = lbegin
c
c  generate the next level by finding all the masked
c  neighbors of nodes in the current level.
c
      do 400 i = lbegin, lvlend
        node = ls(i)
        jstrt = xadj(node)
        jstop = xadj(node + 1) - 1
        if ( jstop .lt. jstrt ) go to 400
        do 300 j = jstrt, jstop
          nbr = adjncy(j)
          if (mask(nbr) .eq. 0) go to 300
          ccsize = ccsize + 1
          ls(ccsize) = nbr
          mask(nbr) = 0
300     continue
400   continue
c
c  compute the current level width.
c  if it is nonzero, generate the next level.
c
      lvsize = ccsize - lvlend
      if (lvsize .gt. 0 ) go to 200
c
c  reset mask to one for the nodes in the level structure.
c
      xls(nlvl+1) = lvlend + 1
      do i = 1, ccsize
        node = ls(i)
        mask(node) = 1
      end do

      return
      end

      
      subroutine fnroot ( root, xadj, adjncy, mask, nlvl, xls, ls)

c*********************************************************************72
c
cc FNROOT implements a modified version of the
c  scheme by gibbs, poole, and stockmeyer to find pseudo-
c  peripheral nodes it determines such a node for the
c  section subgraph specified by mask and root.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters -
c
c  (xadj, adjncy) - adjacency structure pair for the graph.
c  mask - specifies a section subgraph. nodes for which
c  mask is zero are ignored by fnroot.
c
c  updated parameter
c  root - on input, it (along with mask) defines the
c  component for which a pseudo-peripheral node is
c  to be found. on output, it is the node obtained
c
c  output parameters 
c  nlvl - is the number of levels in the level structure
c  rooted at the node root.
c  (xls,ls) - the level structure array pair containing
c  the level structure found.
c
      implicit none
      
      integer adjncy(*)
      integer ls(*)
      integer mask(*), xls(*)
      integer xadj(*), ccsize, j, jstrt, k, kstop, kstrt
      integer mindeg, nabor, ndeg, nlvl, node, nunlvl, root
c
c  determine the level structure rooted at root.
c
      call rootls ( root, xadj, adjncy, mask, nlvl, xls, ls)
      ccsize = xls(nlvl+1) - 1

      if ( nlvl .eq. 1 .or. nlvl .eq. ccsize ) then
        return
      end if
c
c  pick a node with minimum degree from the last level.
c
100   continue

      jstrt = xls(nlvl)
      mindeg = ccsize
      root = ls(jstrt)
      if ( ccsize .eq. jstrt ) go to 400

      do j = jstrt, ccsize

        node = ls(j)
        ndeg = 0
        kstrt = xadj(node)
        kstop = xadj(node+1) - 1

        do k = kstrt, kstop
          nabor = adjncy(k)
          if ( mask(nabor) .gt. 0 ) then
            ndeg = ndeg + 1
          end if
        end do

        if ( ndeg .lt. mindeg ) then
          root = node
          mindeg = ndeg
        end if

      end do
c
c  and generate its rooted level structure.
c
400   continue

      call rootls ( root, xadj, adjncy, mask, nunlvl, xls, ls)
      if (nunlvl .le. nlvl) return
      nlvl = nunlvl
      if (nlvl .lt. ccsize ) go to 100
      return
      end


      subroutine degree ( root, xadj, adjncy, mask, deg, ccsize, ls )

c*********************************************************************72
c
c     DEGREE computes node degrees in a connected component, for the RCM
c     method.
c
c  Discussion:
c
c    This routine computes the degrees of the nodes in the connected
c    component specified by MASK and ROOT.
c
c    Nodes for which MASK is zero are ignored.
c
c  Modified:
c
c    01 January 2009
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  Parameters:
c
c    Input, integer ROOT, the node that defines the component.
c
c    (xadj, adjncy) - adjacency structure pair.
c
c    Input, integer MASK(N), specifies a section subgraph.
c
c    Output, integer DEGREE(N), the degrees of the nodes in the component.
c
c    Output, integer CCSIZE, the size of the component.
c
c    Workspace, integer LS(N), used to store the nodes of the
c    component level by level.
c
      implicit none
      
      integer adjncy(*)
      integer ccsize
      integer deg(*)
      integer i
      integer ideg
      integer j
      integer jstop
      integer jstrt
      integer lbegin
      integer ls(*)
      integer lvlend
      integer lvsize
      integer mask(*)
      integer nbr
      integer node
      integer root
      integer xadj(*)
c
c  Initialization.
c
c  the array xadj is used as a temporary marker to
c  indicate which nodes have been considered so far
c
      ls(1) = root
      xadj(root) = -xadj(root)
      lvlend = 0
      ccsize = 1
c
c  lbegin is the pointer to the beginning of the curi
c  level, and lvlend points to the end of this level.
c
100   continue

      lbegin = lvlend + 1
      lvlend = ccsize
c
c  find the degrees of nodes in the current level,
c  and at the same time, generate the next level.
c
      do i = lbegin, lvlend
        node = ls(i)
        jstrt = -xadj(node)
        jstop = iabs(xadj(node+ 1) ) - 1
        ideg = 0
        if ( jstop .lt. jstrt ) go to 300
        do 200 j = jstrt, jstop
          nbr = adjncy(j)
          if ( mask(nbr) .eq. 0 ) go to 200
          ideg = ideg + 1
          if ( xadj(nbr) .lt. 0 ) go to 200
          xadj(nbr) = -xadj(nbr)
          ccsize = ccsize + 1
          ls(ccsize) = nbr
200     continue
300     deg(node) = ideg
      end do
c
c  Compute the current level width.
c  if it is nonzero , generate another level.
c
      lvsize = ccsize - lvlend
      if ( lvsize .gt. 0 ) go to 100
c
c  Reset xadj to its correct sign and return.
c
      do i = 1, ccsize
        node = ls(i)
        xadj(node) = - xadj(node)
      end do

      return
      end

      
      subroutine rcm ( root, xadj, adjncy, mask, perm, ccsize, deg )

c*********************************************************************72
c
cc RCM renumbers a connected component by the reverse Cuthill McKee
c  algorithm.
c
c  Discussion:
c
c    The connected component is specified by mask and root.
c    The numbering is done starting at the node ROOT, and proceeding
c    using the Reverse Cuthill-McKee algorithm.
c
c    An outline of the algorithm is as follows:
c
c    X(1) = ROOT.
c
c     for ( I = 1 to N-1)
c        Find all unlabeled neighbors of X(I), assign
c        them the next available labels, in order of increasing degree.
c
c    When done, reverse the ordering.
c
c  Modified:
c
c    31 December 2008
c
c  Author:
c
c    Alan George, Joseph Liu
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  input parameters 
c
c    Input, integer ROOT, the node that defines the connected component.
c    It is used as the starting point for the RCM ordering.
c
c    Input, integer XADJ(N+1).  Information about row I is stored
c    in entries XADJ(I) through XADJ(I+1)-1 of ADJNCY.
c
c    Input, integer ADJNCY(*), the adjacency structure. 
c    For each row, the column indices of the nonzero entries.
c
c    Input/output, integer MASK(N), a mask for the nodes.  Only those nodes 
c    with nonzero input mask values are considered by the routine.  The
c    nodes numbered by RCM will have their mask values set to zero.
c
c    Output, integer PERM(N), the RCM ordering.
c
c    Output, integer CCSZE, the size of the connected component
c    that has been numbered.
c
c    Workspace, integer DEG(N), a temporary vector used to hold the degree
c    of the nodes in the section graph specified by MASK and ROOT.
c
      implicit none
      
      integer adjncy(*)
      integer ccsize
      integer deg(*)
      integer fnbr
      integer i
      integer j
      integer jstop
      integer jstrt
      integer k
      integer l
      integer lbegin
      integer lnbr
      integer lperm
      integer lvlend
      integer mask(*)
      integer nbr
      integer node
      integer perm(*)
      integer root
      integer xadj(*)
c
c  Find the degrees of the nodes in the component specified by MASK and ROOT.
c
      call degree ( root, xadj, adjncy, mask, deg, ccsize, perm )

      mask(root) = 0

      if ( ccsize .le. 1 ) then
        return
      end if

      lvlend = 0
      lnbr = 1
c
c  lbegin and lvlend point to the beginning and
c  the end of the current level respectively.
c
100   continue

      lbegin = lvlend + 1
      lvlend = lnbr

      do i = lbegin, lvlend
c
c  For each node in current level.
c
        node = perm(i)
        jstrt = xadj(node)
c       *** NB NEXT LINE WAS IN ERROR ***
c       jstop = xadj(node+i) - 1
        jstop = xadj(node+1) - 1        
c
c  Find the unnumbered neighbors of node.
c
c  fnbr and lnbr point to the first and last
c  unnumbered neighbors respectively of the current
c  node in perm.        
c
        fnbr = lnbr + 1

        do j = jstrt, jstop
          nbr = adjncy(j)
          if ( mask(nbr) .ne. 0 ) then
            lnbr = lnbr + 1
            mask(nbr) = 0
            perm(lnbr) = nbr
          end if
        end do
c
c  sort the neighbors of node in increasing
c  order by degree. linear insertion is used
c
        if ( fnbr .lt. lnbr ) then

          k = fnbr

300       continue

          l = k
          k = k + 1
          nbr = perm(k)

400       if ( l .lt. fnbr ) go to 500
          lperm = perm(l)
          if(deg(lperm) .le. deg(nbr) ) go to 500
          perm(l+1)=lperm
          l=l-1
          go to 400

500       continue

          perm(l+1)=nbr
 
          if ( k .lt. lnbr ) then
            go to 300
          end if

        end if

      end do

      if ( lnbr .gt. lvlend ) then
        go to 100
      end if
c
c  we now have the cuthill mckee ordering.
c  reverse it.
c
      k = ccsize / 2
      l = ccsize

      do i=1, k
        lperm = perm(l)
        perm(l) = perm(i)
        perm(i) = lperm
        l = l - 1
      end do

      return
      end

      
      subroutine i4vec_reverse ( n, a )

c*********************************************************************72
c
cc I4VEC_REVERSE reverses the elements of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    In FORTRAN90, call I4VEC_REVERSE is equivalent to:
c
c      A(1:N) = A(N:1:-1)
c
c  Example:
c
c    Input:
c
c      N = 5,
c      A = ( 11, 12, 13, 14, 15 ).
c
c    Output:
c
c      A = ( 15, 14, 13, 12, 11 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N), the array to be reversed.
c
      implicit none

      integer n
      integer a(n)
      integer i
      integer t

      do i = 1, n / 2
        t        = a(i)
        a(i)     = a(n+1-i)
        a(n+1-i) = t
      end do

      return
      end

      
      subroutine gennd ( neqns, xadj, adjncy, mask,
     &  perm, xls, ls )

c*********************************************************************72
c
cc GENND finds a nested dissection ordering for a general graph.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c    Input, integer NEQNS, the number of equations.
c
c  (xadj, adjncy) - adjacency structure pair.
c
c  output parameters 
c  perm - the nested dissection ordering.
c  
c  working parameters 
c  mask - is used to mask off variables that have
c  been numbered during the orderng process,
c  (xls, ls) - this level structure pair is used as
c  temporary storage by fnroot.
c
      implicit none
      
      integer adjncy(*)
      integer ls(*)
      integer mask(*), perm(*)
      integer xls(*)
      integer xadj(*), i, neqns, nsep, num, root

      do i = 1, neqns
        mask(i) = 1
      end do

      num = 0

      do 300 i = 1, neqns
c
c  for each masked component ...
c
200     if ( mask(i) .eq. 0 ) go to 300
        root = i
c
c  find a separator and number the nodes next,
c
        call fndsep ( root, xadj, adjncy, mask,
     &    nsep, perm(num+1), xls, ls )
        num = num + nsep
        if ( num .ge. neqns ) go to 400
        go to 200
300   continue
c
c  since separators found first should be ordered
c  last, routine revrse is called to adjust the
c  ordering vector,
c
400   continue

      call i4vec_reverse ( neqns, perm )

      return
      end

      
      subroutine fndsep ( root, xadj, adjncy, mask, nsep, sep,
     &  xls, ls )

c*********************************************************************72
c
cc FNDSEP finds a small separator for a connected component.
c
c  Discussion:
c
c    The connected component is specified by MASK.
c
c  Modified:
c
c    24 July 2012 - a number of errors fixed.
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  Parameters:
c
c    Input, integer ROOT, the node that determines the masked component.
c
c    Input, integer XADJ(*), integer ADJNCY(*), the adjacency structure
c    pair.
c
c    Input/output, integer MASK(*), nodes in the separator will have
c    their mask values set to zero.
c
c    Output, integer NSEP, the number of variables in the separator.
c
c    Output, integer SEP(*), the separator nodes.
c
c    Workspace, integer XLS(*), LS(*), the level structure pair for
c    level structure found by FNROOT.
c
      implicit none

      integer adjncy(*)
      integer i
      integer j
      integer jstop
      integer jstrt
      integer ls(*)
      integer mask(*)
      integer midbeg
      integer midend
      integer midlvl
      integer mp1beg
      integer mp1end
      integer nbr
      integer nlvl
      integer node
      integer nsep
      integer root
      integer sep(*)
      integer xls(*)
      integer xadj(*)
c
c  Determine the level structure associated with ROOT.
c
      call rootls ( root, xadj, adjncy, mask, nlvl, xls, ls )
c
c  If the number of levels is less than 3, return the
c  whole component as the separator.
c
      if ( nlvl .lt. 3 ) then

        nsep = xls(nlvl+1) - 1
        do i = 1, nsep
          node = ls(i)
          sep(i) = node
          mask(node) = 0
        end do

        return

      end if
c
c  Find the middle level of the rooted level structure.
c
      midlvl = ( nlvl + 2 ) / 2
      midbeg = xls(midlvl)
      mp1beg = xls(midlvl + 1)
      midend = mp1beg - 1
      mp1end = xls(midlvl+2)
c
c  The separator is obtained by including only those
c  MIDDLE-level nodes with neighbors in the MIDDLE+1
c  level.  XADJ is used temporarily to mark those
c  nodes in the MIDDLE+1 level.
c
      do i = mp1beg, mp1end
        node = ls(i)
        xadj(node) = - xadj(node)
      end do

      nsep = 0
      
      do i = midbeg, midend
        node = ls(i)
        jstrt = xadj(node)
        jstop = iabs ( xadj(node+1) ) - 1

        do j = jstrt, jstop

          nbr = adjncy(j)

          if ( xadj(nbr) .le. 0 ) then
            nsep = nsep + 1
            sep(nsep) = node
            mask(node) = 0
            go to 10
          end if

        end do

10      continue

      end do
c
c  Reset XADJ to its correct sign.
c
      do i = mp1beg, mp1end
        node = ls(i)
        xadj(node) = - xadj(node)
      end do

      return
      end

#else

      
C     ##################################################################
C     3. SPARSPAK.
C     SPARSPAK is no longer (circa 2010) available through netlib.
C      
C     See original descriptions in George & Liu (1981),
C     Computer Solution of Large Sparse Positive-Definite Systems.
C     or (available as PDF) George, Lui & Ng (1994),
C     Computer Solution of Sparse Linear Systems.      
C     ##################################################################
C
C     Here we supply routines from SPARSPAK that:
C     A. Carry out Reverse Cuthill--McKee node ordering (Ch 4 of book).
C        Main routine is RCM
C     B. Carry out nested dissection ordering           (Ch 8 of book).
C        Main routine is GENND
C
C     Versions below in UPPERCASE are as per book (with fixes).

C----- SUBROUTINE ROOTLS
C***************************************************************        
C***************************************************************        
C********     ROOTLS ..... ROOTED LEVEL STRUCTURE      *********        
C***************************************************************        
C***************************************************************        
C                                                                       
C     PURPOSE - ROOTLS GENERATES THE LEVEL STRUCTURE ROOTED             
C        AT THE INPUT NODE CALLED ROOT. ONLY THOSE NODES FOR            
C        WHICH MASK IS NONZERO WILL BE CONSIDERED.                      
C                                                                       
C     INPUT PARAMETERS -                                                
C        ROOT - THE NODE AT WHICH THE LEVEL STRUCTURE IS TO             
C               BE ROOTED.                                              
C        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR FOR THE              
C               GIVEN GRAPH.                                            
C        MASK - IS USED TO SPECIFY A SECTION SUBGRAPH. NODES            
C               WITH MASK(I)=0 ARE IGNORED.                             
C                                                                       
C     OUTPUT PARAMETERS -                                               
C        NLVL - IS THE NUMBER OF LEVELS IN THE LEVEL STRUCTURE.         
C        (XLS, LS) - ARRAY PAIR FOR THE ROOTED LEVEL STRUCTURE.         
C                                                                       
C***************************************************************        
C                                                                       
      SUBROUTINE  ROOTLS ( ROOT, XADJ, ADJNCY, MASK, NLVL, XLS, LS )    
C                                                                       
C***************************************************************        
C
      IMPLICIT NONE
      
      INTEGER ADJNCY(1), LS(1), MASK(1), XLS(1)                      
      INTEGER XADJ(1), I, J, JSTOP, JSTRT, LBEGIN,                   
     1     CCSIZE, LVLEND, LVSIZE, NBR, NLVL,                     
     1     NODE, ROOT                                             
C                                                                       
C***************************************************************        
C                                                                       
C        ------------------                                             
C        INITIALIZATION ...                                             
C        ------------------                                             
         MASK(ROOT) = 0                                                 
         LS(1) = ROOT                                                   
         NLVL = 0                                                       
         LVLEND = 0                                                     
         CCSIZE = 1                                                     
C        -----------------------------------------------------          
C        LBEGIN IS THE POINTER TO THE BEGINNING OF THE CURRENT          
C        LEVEL, AND LVLEND POINTS TO THE END OF THIS LEVEL.             
C        -----------------------------------------------------          
  200    LBEGIN = LVLEND + 1                                            
         LVLEND = CCSIZE                                                
         NLVL = NLVL + 1                                                
         XLS(NLVL) = LBEGIN                                             
C        -------------------------------------------------              
C        GENERATE THE NEXT LEVEL BY FINDING ALL THE MASKED              
C        NEIGHBORS OF NODES IN THE CURRENT LEVEL.                       
C        -------------------------------------------------              
         DO 400 I = LBEGIN, LVLEND                                      
            NODE = LS(I)                                                
            JSTRT = XADJ(NODE)                                          
            JSTOP = XADJ(NODE + 1) - 1                                  
            IF ( JSTOP .LT. JSTRT )  GO TO 400                          
               DO 300 J = JSTRT, JSTOP                                  
                  NBR = ADJNCY(J)                                       
                  IF (MASK(NBR) .EQ. 0) GO TO 300                       
                     CCSIZE = CCSIZE + 1                                
                     LS(CCSIZE) = NBR                                   
                     MASK(NBR) = 0                                      
  300          CONTINUE                                                 
  400    CONTINUE                                                       
C        ------------------------------------------                     
C        COMPUTE THE CURRENT LEVEL WIDTH.                               
C        IF IT IS NONZERO, GENERATE THE NEXT LEVEL.                     
C        ------------------------------------------                     
         LVSIZE = CCSIZE - LVLEND                                       
         IF (LVSIZE .GT. 0 ) GO TO 200                                  
C        -------------------------------------------------------        
C        RESET MASK TO ONE FOR THE NODES IN THE LEVEL STRUCTURE.        
C        -------------------------------------------------------        
         XLS(NLVL+1) = LVLEND + 1                                       
         DO 500 I = 1, CCSIZE                                           
            NODE = LS(I)                                                
            MASK(NODE) = 1                                              
  500    CONTINUE                                                       
         RETURN                                                         
      END                                                               
                                                                        
                                                                        
C----- SUBROUTINE FNROOT                                                
C***************************************************************        
C***************************************************************        
C*******     FNROOT ..... FIND PSEUDO-PERIPHERAL NODE    *******        
C***************************************************************        
C***************************************************************        
C                                                                       
C    PURPOSE - FNROOT IMPLEMENTS A MODIFIED VERSION OF THE              
C       SCHEME BY GIBBS, POOLE, AND STOCKMEYER TO FIND PSEUDO-          
C       PERIPHERAL NODES.  IT DETERMINES SUCH A NODE FOR THE            
C       SECTION SUBGRAPH SPECIFIED BY MASK AND ROOT.                    
C                                                                       
C    INPUT PARAMETERS -                                                 
C       (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR FOR THE GRAPH.        
C       MASK - SPECIFIES A SECTION SUBGRAPH. NODES FOR WHICH            
C              MASK IS ZERO ARE IGNORED BY FNROOT.                      
C                                                                       
C    UPDATED PARAMETER -                                                
C       ROOT - ON INPUT, IT (ALONG WITH MASK) DEFINES THE               
C              COMPONENT FOR WHICH A PSEUDO-PERIPHERAL NODE IS          
C              TO BE FOUND. ON OUTPUT, IT IS THE NODE OBTAINED.         
C                                                                       
C    OUTPUT PARAMETERS -                                                
C       NLVL - IS THE NUMBER OF LEVELS IN THE LEVEL STRUCTURE           
C              ROOTED AT THE NODE ROOT.                                 
C       (XLS,LS) - THE LEVEL STRUCTURE ARRAY PAIR CONTAINING            
C                  THE LEVEL STRUCTURE FOUND.                           
C                                                                       
C    PROGRAM SUBROUTINES -                                              
C       ROOTLS.                                                         
C                                                                       
C***************************************************************        
C                                                                       
      SUBROUTINE  FNROOT ( ROOT, XADJ, ADJNCY, MASK, NLVL, XLS, LS )    
C                                                                       
C***************************************************************        
C
      IMPLICIT NONE
      
         INTEGER ADJNCY(1), LS(1), MASK(1), XLS(1)                      
         INTEGER XADJ(1), CCSIZE, J, JSTRT, K, KSTOP, KSTRT,            
     1           MINDEG, NABOR, NDEG, NLVL, NODE, NUNLVL,               
     1           ROOT                                                   
C                                                                       
C***************************************************************        
C                                                                       
C        ---------------------------------------------                  
C        DETERMINE THE LEVEL STRUCTURE ROOTED AT ROOT.                  
C        ---------------------------------------------                  
         CALL  ROOTLS ( ROOT, XADJ, ADJNCY, MASK, NLVL, XLS, LS )       
         CCSIZE = XLS(NLVL+1) - 1                                       
         IF ( NLVL .EQ. 1 .OR. NLVL .EQ. CCSIZE ) RETURN                
C        ----------------------------------------------------           
C        PICK A NODE WITH MINIMUM DEGREE FROM THE LAST LEVEL.           
C        ----------------------------------------------------           
  100    JSTRT = XLS(NLVL)                                              
         MINDEG = CCSIZE                                                
         ROOT = LS(JSTRT)                                               
         IF ( CCSIZE .EQ. JSTRT )  GO TO 400                            
            DO 300 J = JSTRT, CCSIZE                                    
               NODE = LS(J)                                             
               NDEG = 0                                                 
               KSTRT = XADJ(NODE)                                       
               KSTOP = XADJ(NODE+1) - 1                                 
               DO 200 K = KSTRT, KSTOP                                  
                  NABOR = ADJNCY(K)                                     
                  IF ( MASK(NABOR) .GT. 0 )  NDEG = NDEG + 1            
  200          CONTINUE                                                 
               IF ( NDEG .GE. MINDEG ) GO TO 300                        
                  ROOT = NODE                                           
                  MINDEG = NDEG                                         
  300       CONTINUE                                                    
C        ----------------------------------------                       
C        AND GENERATE ITS ROOTED LEVEL STRUCTURE.                       
C        ----------------------------------------                       
  400    CALL  ROOTLS ( ROOT, XADJ, ADJNCY, MASK, NUNLVL, XLS, LS )     
         IF (NUNLVL .LE. NLVL)  RETURN                                  
            NLVL = NUNLVL                                               
            IF ( NLVL .LT. CCSIZE )  GO TO 100                          
            RETURN                                                      
      END                                                               

                                                                        
C----- SUBROUTINE DEGREE                                                
C***************************************************************        
C***************************************************************        
C********     DEGREE ..... DEGREE IN MASKED COMPONENT   ********        
C***************************************************************        
C***************************************************************        
C                                                                       
C     PURPOSE - THIS ROUTINE COMPUTES THE DEGREES OF THE NODES          
C        IN THE CONNECTED COMPONENT SPECIFIED BY MASK AND ROOT.         
C        NODES FOR WHICH MASK IS ZERO ARE IGNORED.                      
C                                                                       
C     INPUT PARAMETER -                                                 
C        ROOT - IS THE INPUT NODE THAT DEFINES THE COMPONENT.           
C        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR.                     
C        MASK - SPECIFIES A SECTION SUBGRAPH.                           
C                                                                       
C     OUTPUT PARAMETERS -                                               
C        DEG - ARRAY CONTAINING THE DEGREES OF THE NODES IN             
C              THE COMPONENT.                                           
C        CCSIZE-SIZE OF THE COMPONENT SPECIFED BY MASK AND ROOT         
C                                                                       
C     WORKING PARAMETER -                                               
C        LS - A TEMPORARY VECTOR USED TO STORE THE NODES OF THE         
C               COMPONENT LEVEL BY LEVEL.                               
C                                                                       
C***************************************************************        
C                                                                       
      SUBROUTINE  DEGREE ( ROOT, XADJ, ADJNCY, MASK,                    
     1                     DEG, CCSIZE, LS )                            
C                                                                       
C***************************************************************        
C
      IMPLICIT NONE
      
         INTEGER ADJNCY(1), DEG(1), LS(1), MASK(1)                      
         INTEGER XADJ(1), CCSIZE, I, IDEG, J, JSTOP, JSTRT,             
     1           LBEGIN, LVLEND, LVSIZE, NBR, NODE, ROOT                
C                                                                       
C***************************************************************        
C                                                                       
C        -------------------------------------------------              
C        INITIALIZATION ...                                             
C        THE ARRAY XADJ IS USED AS A TEMPORARY MARKER TO                
C        INDICATE WHICH NODES HAVE BEEN CONSIDERED SO FAR.              
C        -------------------------------------------------              
         LS(1) = ROOT                                                   
         XADJ(ROOT) = -XADJ(ROOT)                                       
         LVLEND = 0                                                     
         CCSIZE = 1                                                     
C        -----------------------------------------------------          
C        LBEGIN IS THE POINTER TO THE BEGINNING OF THE CURRENT          
C        LEVEL, AND LVLEND POINTS TO THE END OF THIS LEVEL.             
C        -----------------------------------------------------          
  100    LBEGIN = LVLEND + 1                                            
         LVLEND = CCSIZE                                                
C        -----------------------------------------------                
C        FIND THE DEGREES OF NODES IN THE CURRENT LEVEL,                
C        AND AT THE SAME TIME, GENERATE THE NEXT LEVEL.                 
C        -----------------------------------------------                
         DO 400 I = LBEGIN, LVLEND                                      
            NODE = LS(I)                                                
            JSTRT = -XADJ(NODE)                                         
            JSTOP = IABS(XADJ(NODE + 1)) - 1                            
            IDEG = 0                                                    
            IF ( JSTOP .LT. JSTRT ) GO TO 300                           
               DO 200 J = JSTRT, JSTOP                                  
                  NBR = ADJNCY(J)                                       
                  IF ( MASK(NBR) .EQ. 0 )  GO TO  200                   
                     IDEG = IDEG + 1                                    
                     IF ( XADJ(NBR) .LT. 0 ) GO TO 200                  
                        XADJ(NBR) = -XADJ(NBR)                          
                        CCSIZE = CCSIZE + 1                             
                        LS(CCSIZE) = NBR                                
  200          CONTINUE                                                 
  300       DEG(NODE) = IDEG                                            
  400    CONTINUE                                                       
C        ------------------------------------------                     
C        COMPUTE THE CURRENT LEVEL WIDTH.                               
C        IF IT IS NONZERO , GENERATE ANOTHER LEVEL.                     
C        ------------------------------------------                     
         LVSIZE = CCSIZE - LVLEND                                       
         IF ( LVSIZE .GT. 0 ) GO TO 100                                 
C        ------------------------------------------                     
C        RESET XADJ TO ITS CORRECT SIGN AND RETURN.                     
C        ------------------------------------------                     
         DO 500 I = 1, CCSIZE                                           
            NODE = LS(I)                                                
            XADJ(NODE) = -XADJ(NODE)                                    
  500    CONTINUE                                                       
         RETURN                                                         
      END                                                               

                                                                        
C----- SUBROUTINE RCM                                                   
C***************************************************************        
C***************************************************************        
C********     RCM ..... REVERSE CUTHILL-MCKEE ORDERING   *******        
C***************************************************************        
C***************************************************************        
C                                                                       
C     PURPOSE - RCM NUMBERS A CONNECTED COMPONENT SPECIFIED BY          
C        MASK AND ROOT, USING THE RCM ALGORITHM.                        
C        THE NUMBERING IS TO BE STARTED AT THE NODE ROOT.               
C                                                                       
C     INPUT PARAMETERS -                                                
C        ROOT - IS THE NODE THAT DEFINES THE CONNECTED                  
C               COMPONENT AND IT IS USED AS THE STARTING                
C               NODE FOR THE RCM ORDERING.                              
C        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR FOR                  
C               THE GRAPH.                                              
C                                                                       
C     UPDATED PARAMETERS -                                              
C        MASK - ONLY THOSE NODES WITH NONZERO INPUT MASK                
C               VALUES ARE CONSIDERED BY THE ROUTINE.  THE              
C               NODES NUMBERED BY RCM WILL HAVE THEIR                   
C               MASK VALUES SET TO ZERO.                                
C                                                                       
C     OUTPUT PARAMETERS -                                               
C        PERM - WILL CONTAIN THE RCM ORDERING.                          
C        CCSIZE - IS THE SIZE OF THE CONNECTED COMPONENT                
C               THAT HAS BEEN NUMBERED BY RCM.                          
C                                                                       
C     WORKING PARAMETER -                                               
C        DEG - IS A TEMPORARY VECTOR USED TO HOLD THE DEGREE            
C               OF THE NODES IN THE SECTION GRAPH SPECIFIED             
C               BY MASK AND ROOT.                                       
C                                                                       
C     PROGRAM SUBROUTINES -                                             
C        DEGREE.                                                        
C                                                                       
C***************************************************************        
C                                                                       
      SUBROUTINE  RCM ( ROOT, XADJ, ADJNCY, MASK,                       
     1                  PERM, CCSIZE, DEG )                             
C                                                                       
C***************************************************************        
C
      IMPLICIT NONE
      
         INTEGER ADJNCY(1), DEG(1), MASK(1), PERM(1)                    
         INTEGER XADJ(1), CCSIZE, FNBR, I, J, JSTOP,                    
     1           JSTRT, K, L, LBEGIN, LNBR, LPERM,                      
     1           LVLEND, NBR, NODE, ROOT                                
C                                                                       
C***************************************************************        
C                                                                       
C        -------------------------------------                          
C        FIND THE DEGREES OF THE NODES IN THE                           
C        COMPONENT SPECIFIED BY MASK AND ROOT.                          
C        -------------------------------------                          
         CALL  DEGREE ( ROOT, XADJ, ADJNCY, MASK, DEG,                  
     1                  CCSIZE, PERM )                                  
         MASK(ROOT) = 0                                                 
         IF ( CCSIZE .LE. 1 ) RETURN                                    
         LVLEND = 0                                                     
         LNBR = 1                                                       
C        --------------------------------------------                   
C        LBEGIN AND LVLEND POINT TO THE BEGINNING AND                   
C        THE END OF THE CURRENT LEVEL RESPECTIVELY.                     
C        --------------------------------------------                   
  100    LBEGIN = LVLEND + 1                                            
         LVLEND = LNBR                                                  
         DO 600 I = LBEGIN, LVLEND                                      
C           ----------------------------------                          
C           FOR EACH NODE IN CURRENT LEVEL ...                          
C           ----------------------------------                          
            NODE = PERM(I)                                              
            JSTRT = XADJ(NODE)                                          
            JSTOP = XADJ(NODE+1) - 1                                    
C           ------------------------------------------------            
C           FIND THE UNNUMBERED NEIGHBORS OF NODE.                      
C           FNBR AND LNBR POINT TO THE FIRST AND LAST                   
C           UNNUMBERED NEIGHBORS RESPECTIVELY OF THE CURRENT            
C           NODE IN PERM.                                               
C           ------------------------------------------------            
            FNBR = LNBR + 1                                             
            DO 200 J = JSTRT, JSTOP                                     
               NBR = ADJNCY(J)                                          
               IF ( MASK(NBR) .EQ. 0 )  GO TO 200                       
                  LNBR = LNBR + 1                                       
                  MASK(NBR) = 0                                         
                  PERM(LNBR) = NBR                                      
  200       CONTINUE                                                    
            IF ( FNBR .GE. LNBR )  GO TO 600                            
C              ------------------------------------------               
C              SORT THE NEIGHBORS OF NODE IN INCREASING                 
C              ORDER BY DEGREE. LINEAR INSERTION IS USED.               
C              ------------------------------------------               
               K = FNBR                                                 
  300          L = K                                                    
                  K = K + 1                                             
                  NBR = PERM(K)                                         
  400             IF ( L .LT. FNBR )  GO TO 500                         
                     LPERM = PERM(L)                                    
                     IF ( DEG(LPERM) .LE. DEG(NBR) )  GO TO 500         
                        PERM(L+1) = LPERM                               
                        L = L - 1                                       
                        GO TO 400                                       
  500             PERM(L+1) = NBR                                       
                  IF ( K .LT. LNBR )  GO TO 300                         
  600    CONTINUE                                                       
         IF (LNBR .GT. LVLEND) GO TO 100                                
C        ---------------------------------------                        
C        WE NOW HAVE THE CUTHILL MCKEE ORDERING.                        
C        REVERSE IT BELOW ...                                           
C        ---------------------------------------                        
         K = CCSIZE/2                                                   
         L = CCSIZE                                                     
         DO 700 I = 1, K                                                
            LPERM = PERM(L)                                             
            PERM(L) = PERM(I)                                           
            PERM(I) = LPERM                                             
            L = L - 1                                                   
  700    CONTINUE                                                       
         RETURN                                                         
      END                                                               


C----- SUBROUTINE REVRSE
C***************************************************************
C***************************************************************      
C*********    REVRSE .... REVERSE ORDERING of VECTOR   *********
C***************************************************************
C***************************************************************
C
C     THIS IS A UTILITY ROUTINE DESCRIBED BUT NOT SUPPLIED IN
C     SECTION 7.3.2 OF GEORGE & LIU
C
C***************************************************************
C      
      SUBROUTINE REVRSE ( NV, V )
C
      IMPLICIT NONE
      
      INTEGER V(1), I, NV, NV2, TMP
C
      NV2 = NV / 2
C
      DO 100 I = 1, NV2
         TMP           = V(I)
         V(I)          = V(NV - I + 1)
         V(NV - I + 1) = TMP
 100  CONTINUE
      RETURN
      END
      
C----- SUBROUTINE GENND
C***************************************************************
C***************************************************************
C*********    GENND ..... GENERAL NESTED DISSECTION     ********
C***************************************************************
C***************************************************************
C     
C     PURPOSE - SUBROUTINE GENND FINDS A NESTED DISSECTION
C     ORDERING FOR A GENERAL GRAPH.
C     
C     
C     INPUT PARAMETERS -
C     NEQNS - NUMBER OF EQUATIONS.
C     (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR.
C     
C     OUTPUT PARAMETERS -
C     PERM - THE NESTED DISSECTION ORDERING.
C     
C     WORKING PARAMETERS -
C     MASK - IS USED TO MASK OFF VARIABLES THAT HAVE
C     BEEN NUMBERED DURING THE ORDERNG PROCESS.
C     (XLS, LS) - THIS LEVEL STRUCTURE PAIR IS USED AS
C     TEMPORARY STORAGE BY FNROOT.
C     
C     PROGRAM SUBROUTINES -
C     FNDSEP, REVRSE.
C     
C***************************************************************
C     
      SUBROUTINE GENND ( NEQNS, XADJ, ADJNCY, MASK,
     1                   PERM, XLS, LS )
C     
C***************************************************************
C
      IMPLICIT NONE
      
      INTEGER ADJNCY(1), MASK(1), LS(1), PERM(1),
     1        XLS(1)
      INTEGER XADJ(1), I, NEQNS, NSEP, NUM, ROOT
C     
C***************************************************************
C

      DO 100 I = 1, NEQNS
         MASK(I) = 1
 100  CONTINUE
      NUM = 0
      DO 300 I = 1, NEQNS
C        -----------------------------
C        FOR EACH MASKED COMPONENT ...
C        -----------------------------
 200     IF ( MASK(I) .EQ. 0 ) GO TO 300
         ROOT = I
C        -------------------------------------------
C        FIND A SEPARATOR AND NUMBER THE NODES NEXT.
C     -------------------------------------------

         CALL FNDSEP ( ROOT, XADJ, ADJNCY, MASK,
     1     NSEP, PERM(NUM+1), XLS, LS )
           NUM  = NUM + NSEP

         IF ( NUM .GE. NEQNS ) GO TO 400

         GO TO 200
 300  CONTINUE
C     ----------------------------------------------
C     SINCE SEPARATORS FOUND FIRST SHOULD BE ORDERED
C     LAST, ROUTINE REVRSE IS CALLED TO ADJUST THE
C     ORDERING VECTOR
C     ----------------------------------------------
 400  CALL REVRSE ( NEQNS, PERM )
      RETURN
      END


C----- SUBROUTINE FNDSEP      
C***************************************************************
C***************************************************************
C************      FNDSEP ..... FIND SEPARATOR      ************
C***************************************************************
C***************************************************************
C     
C     PURPOSE - THIS ROUTINE IS USED TO FIND A SMALL
C     SEPARATOR FOR A CONNECTED COMPONENT SPECIFIED
C     BY MASK IN THE GIVEN GRAPH.
C     
C     INPUT PARAMETERS -
C     ROOT - IS THE NODE THAT DETERMINES THE MASKED
C     COMPONENT.
C     (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE PAIR.
C     
C     OUTPUT PARAMETERS -
C     NSEP - NUMBER OF VARIABLES IN THE SEPARATOR.
C     SEP - VECTOR CONTAINING THE SEPARATOR NODES.
C     
C     UPDATED PARAMETER -
C     MASK - NODES IN THE SEPARATOR HAVE THEIR MASK
C     VALUES SET TO ZERO.
C     
C     WORKING PARAMETERS -
C     (XLS, LS) - LEVEL STRUCTURE PAIR FOR LEVEL STRUCTURE
C     FOUND BY FNROOT.
C     
C     PROGRAM SUBROUTINES -
C     FNROOT.
C     
C***************************************************************
C     
      SUBROUTINE FNDSEP ( ROOT, XADJ, ADJNCY, MASK,
     1                    NSEP, SEP, XLS , LS )
C     
C***************************************************************
C
      IMPLICIT NONE
      
      INTEGER ADJNCY(1), LS(1), MASK(1), SEP(1), XLS(1)
      INTEGER XADJ(1), I, J, JSTOP, JSTRT, MIDBEG,
     1        MIDEND, MIDLVL, MP1BEG, MP1END,
     1        NBR, NLVL, NODE, NSEP, ROOT
C     
C***************************************************************
C     
C      CALL FNROOT ( ROOT, XADJ, ADJNCY, MASK,
C     1     NLVL, XLS, LS )
      CALL ROOTLS ( ROOT, XADJ, ADJNCY, MASK,
     1              NLVL, XLS, LS )      
C     ----------------------------------------------
C     IF THE NUMBER OF LEVELS IS LESS THAN 3, RETURN
C     THE WHOLE COMPONENT AS THE SEPARATOR.
C     ----------------------------------------------

C     GOTO 200

      IF ( NLVL .GE. 3 ) GO TO 200
      NSEP = XLS(NLVL+1) - 1
      DO 100 I = 1, NSEP
         NODE = LS(I)
         SEP(I) = NODE
         MASK(NODE) = 0
 100  CONTINUE
      RETURN

C     ----------------------------------------------------
C     FIND THE MIDDLE LEVEL OF THE ROOTED LEVEL STRUCTURE.
C     ----------------------------------------------------
 200  MIDLVL = (NLVL + 2)/2
      MIDBEG = XLS(MIDLVL)
C     MPIBEG = XLS(MIDLVL + 1)
      MP1BEG = XLS(MIDLVL + 1)       
      MIDEND = MP1BEG - 1
C     MPIEND = XLS(MIDLVL+2) - 1
      MP1END = XLS(MIDLVL+2)
C     -------------------------------------------------
C     THE SEPARATOR IS OBTAINED BY INCLUDING ONLY THOSE
C     MIDDLE-LEVEL NODES WITH NEIGHBORS IN THE MIDDLE+1
C     LEVEL. XADJ IS USED TEMPORARILY TO MARK THOSE
C     NODES IN THE MIDDLE+1 LEVEL.
C     -------------------------------------------------
C     DO 300 I = MPIBEG, MP1END
      DO 300 I = MP1BEG, MP1END         
         NODE = LS(I)
         XADJ(NODE) = - XADJ(NODE)
 300  CONTINUE
      NSEP  = 0
      DO 500 I = MIDBEG, MIDEND
         NODE = LS(I)
         JSTRT = XADJ(NODE)
         JSTOP = IABS(XADJ(NODE+1)) - 1
         DO 400 J = JSTRT, JSTOP
            NBR = ADJNCY(J)
            IF ( XADJ(NBR) .GT. 0 ) GO TO 400
            NSEP = NSEP + 1
            SEP(NSEP) = NODE
            MASK(NODE) = 0
            GO TO 500
 400     CONTINUE
 500  CONTINUE
C     -------------------------------
C     RESET XADJ TO ITS CORRECT SIGN.
C     -------------------------------
      DO 600 I = MP1BEG, MP1END
         NODE = LS(I)
         XADJ(NODE) = - XADJ(NODE)
 600  CONTINUE
      RETURN
      END

#endif

