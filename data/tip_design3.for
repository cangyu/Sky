$debug
c=======================================================c
c===============  PROGRAM MAIN  ========================c
c=======================================================c

	dimension rootx(400)
	dimension up_y(400,1000),down_y(400,1000)
	integer   number_section(400)

	dimension sect_up_y(400,1000),sect_down_y(400,1000)

	dimension xx(400,1000),yy(400,1000),zz(400,1000)
	integer   num(400)

      common /n_vertex/number_vertex
      common /n_edge  /number_edge
      common /n_face  /number_face

c...   
      number_vertex = 0
	number_edge   = 0 
	number_face   = 0
	nfile         = 6 

	open(nfile,file='point.jou')

ccc.. centre body design:
c...1
	open(3,file='table11.dat')

	number_section(1)= 103
	do i=1,number_section(1)
	read(3,*) rootx(i), sect_up_y(1,i),sect_down_y(1,i)
      
	sect_up_y(1,i)=1.5*sect_up_y(1,i)
	sect_down_y(1,i)=1.5*sect_down_y(1,i)

	end do
	close(3)


c...2
	open(3,file='table13.dat')

	number_section(2)= 103
	do i=1,number_section(2)
	read(3,*) rootx(i), sect_up_y(2,i),sect_down_y(2,i)

	end do
	close(3)
c...3
	open(3,file='table12.dat')
	number_section(3)= 103
	do i=1,number_section(3)
	read(3,*) rootx(i), sect_up_y(3,i),sect_down_y(3,i)
	end do
	close(3)


c...4
	open(3,file='table10.dat')
	number_section(4)= 103
	do i=1,number_section(4)
	read(3,*) rootx(i), sect_up_y(4,i),sect_down_y(4,i)
	end do
	close(3)

c===================================================c
c=======  Design parameter  ========================c


ccc... inner wing design: 1:

c..  design parameter: eb_body is the length, el is spanwise length.
	el_in1 = 650.0
      el_all = 4270.0
	eb_in1 = 2450.0
	ez0    = 0.0
c..  leading edge: maybe the angle is changable.
	ang_le  =  tan(3.1415926*75./180.)

c.. orignal coordinates for central body:leading edge
	x_in1_le  =  0.0
	z_in1_le  =  0.0

c..  trailing edge: maybe the angle is changable.
	ang_te  =  tan(-3.1415926*89./180.)

c..  orignal coordinates for central body: trailing edge
	x_in1_te  =  eb_in1
	z_in1_te  =  0.0

c..  twist angle: twist around leading edge
	ang_twist_w1  =  -2.0
      ang_sf		  =  0.0*3.1415926/180.

c.. design section numbers:
	mmb = 10
	eel = el_in1/real(mmb)

c1
	do i = 1,mmb + 1
        num(i) = number_section(1)
	end do

c... linear interpolate for media airfoil:
	do i =1, number_section(1)
	   ey_up   = sect_up_y  (1,i) - sect_up_y  (2,i)
	   ey_down = sect_down_y(1,i) - sect_down_y(2,i)

	   elup    = ey_up  /real(mmb)
	   eldown  = ey_down/real(mmb)

	   do k =0, mmb
	      up_y  (k+1,i) = sect_up_y  (1,i)  - elup  *real(k)
	      down_y(k+1,i) = sect_down_y(1,i)  - eldown*real(k)
	   end do

	end do

      ex_le0 = x_in1_le
      ex_te0 = x_in1_te

	do k = 0, mmb
	  ez = real(k)*eel 

c... hou lue:  	  
	  zp = ez + z_in1_le
  	  ang_l = 90.-(15.+(650.-zp)/650.*50.)
	  ang_t = -(89.-30.*(650.-zp)/650.)
	  ang_le  =    tan((3.1415926*ang_l)/180.)
  	  ang_te  =    tan((3.1415926*ang_t)/180.)

	  ex_le = ex_le0	+ eel/(ang_le)
	  ex_te	= ex_te0	+ eel/(ang_te)
	  eell	= abs(ex_te - ex_le)

        ex_le0 = ex_le
        ex_te0 = ex_te

       print*, ang_l, ang_t,eell,'mumumu'
	pause

c... twist:
	zp = ez + z_in1_le
	em = (zp-0.0)/(el_all-0.0)
        pai = 3.1415926/180.
	ang_twist = (twist_wing(em)*ang_twist_w1)*pai

c... shang fan:
	  d_y = (zp-0.0)*ang_sf

	  do i=1,num(k+1)

	     xm	 = 	rootx(i)*eell
	     ym	 =  up_y(k+1,i)*eell

	     call twist(xm,ym,ang_twist,xe,ye)

	     xx(k+1,i) = xe + ex_le
	     yy(k+1,i) = ye + d_y
	     zz(k+1,i) = ez + z_in1_le
	  end do

	end do

	call face_create_skin(mmb+1,num,xx,yy,zz,nfile)

      ex_le0 = x_in1_le
      ex_te0 = x_in1_te

	do k = 0, mmb
	  ez = real(k)*eel

c... hou lue:
	  zp = ez + z_in1_le
	ang_l = 90.-(15.+(650.-zp)/650.*50.)
	ang_t = -(89.-30.*(650.-zp)/650.)
	  ang_le  =    tan((3.1415926*ang_l)/180.)
  	  ang_te  =    tan((3.1415926*ang_t)/180.)


	  ex_le = ex_le0 + eel/(ang_le)
	  ex_te	= ex_te0 + eel/(ang_te)
	  eell	= abs(ex_te - ex_le)
	ex_le0 =ex_le
	ex_te0 =ex_te

c... twist:
	zp = ez + z_in1_le
	em = (zp-0.0)/(el_all-0.0)
        pai = 3.1415926/180.
	ang_twist = (twist_wing(em)*ang_twist_w1)*pai

c	ang_twist = twist_wing(em)*3.1415926/180.*ang_twist_w1

c... shang fan:
	  d_y = (zp-0.0)*ang_sf

	  do i=1,num(k+1)

	     xm	 = 	rootx(i)*eell
	     ym	 =  down_y(k+1,i)*eell

	     call twist(xm,ym,ang_twist,xe,ye)

	     xx(k+1,i) = xe + ex_le
	     yy(k+1,i) = ye	+ d_y
	     zz(k+1,i) = ez + z_in1_le
	  end do

	end do

	call face_create_skin(mmb+1,num,xx,yy,zz,nfile)
c	goto 111

ccc... inner wing design: 1:

c..  design parameter: eb_body is the length, el is spanwise length.
	el_in1 = 2720.
	ez0    = ez

c..  leading edge: maybe the angle is changable.
	ang_le  =  tan(3.1415926*75./180.)

c.. orignal coordinates for central body:leading edge
	x_in1_le  =  ex_le
	z_in1_le  =  zp

c..  trailing edge: maybe the angle is changable.
	ang_te  =  tan(-3.1415926*89.0/180.)


c..  orignal coordinates for central body: trailing edge
	x_in1_te  =  ex_te
	z_in1_te  =  zp

c..  twist angle: twist around leading edge
	ang_twist_w1  =  -2.0

      ang_sf		  =  4.0*3.1415926/180.

c.. design section numbers:
	mmb = 10
	eel = el_in1/real(mmb)

c1
	do i = 1,mmb + 1
        num(i) = number_section(1)
	end do

c... linear interpolate for media airfoil:
	do i =1, number_section(1)
	   ey_up   = sect_up_y  (2,i) - sect_up_y  (3,i)
	   ey_down = sect_down_y(2,i) - sect_down_y(3,i)

	   elup    = ey_up  /real(mmb)
	   eldown  = ey_down/real(mmb)

	   do k =0, mmb
	      up_y  (k+1,i) = sect_up_y  (2,i)  - elup  *real(k)
	      down_y(k+1,i) = sect_down_y(2,i)  - eldown*real(k)
	   end do

	end do

	do k = 0, mmb
	  ez = real(k)*eel 
	  ex_le = x_in1_le	+ ez/(ang_le)
	  ex_te	= x_in1_te	+ ez/(ang_te)
	  eell	= abs(ex_te - ex_le)

c... twist:
	zp = ez + z_in1_le

	em = (zp-0.0)/(el_all-0.0)
        pai = 3.1415926/180.
	ang_twist = (twist_wing(em)*ang_twist_w1)*pai

c	ang_twist = twist_wing(em)*3.1415926/180.*ang_twist_w1
        
c        ang_twist = 0.0    
c        ang_sf=0.0
c... shang fan:
	  d_y = (zp-650.0)*ang_sf

	  do i=1,num(k+1)

	     xm	 = 	rootx(i)*eell
	     ym	 =  up_y(k+1,i)*eell

	     call twist(xm,ym,ang_twist,xe,ye)

	     xx(k+1,i) = xe + ex_le
	     yy(k+1,i) = ye + d_y
	     zz(k+1,i) = ez + z_in1_le
	  end do

	end do

	call face_create_skin(mmb+1,num,xx,yy,zz,nfile)

	do k = 0, mmb
	  ez = real(k)*eel
	  ex_le = x_in1_le + ez/(ang_le)
	  ex_te	= x_in1_te + ez/(ang_te)
	  eell	= abs(ex_te - ex_le)

c... twist:
	zp = ez + z_in1_le
	em = (zp-0.0)/(el_all-0.0)
        pai = 3.1415926/180.
	ang_twist = (twist_wing(em)*ang_twist_w1)*pai

c	ang_twist = twist_wing(em)*3.1415926/180.*ang_twist_w1
c        ang_twist = 0.0    
c        ang_sf=0.0
        
c... shang fan:
	  d_y = (zp-650.0)*ang_sf

	  do i=1,num(k+1)

	     xm	 = 	rootx(i)*eell
	     ym	 =  down_y(k+1,i)*eell

	     call twist(xm,ym,ang_twist,xe,ye)

	     xx(k+1,i) = xe + ex_le
	     yy(k+1,i) = ye	+ d_y
	     zz(k+1,i) = ez + z_in1_le
	  end do

	end do

	call face_create_skin(mmb+1,num,xx,yy,zz,nfile)

ccc... inner wing design: 2:

c..  design parameter: eb_body is the length, el is spanwise length.
	el_in1 = 900.
	ez0    = zp

c..  leading edge: maybe the angle is changable.
	ang_le  =  tan(3.1415926*75./180.)

c.. orignal coordinates for central body:leading edge
	x_in1_le  =  ex_le
	z_in1_le  =  zp



c..  trailing edge: maybe the angle is changable.
	ang_te  =  tan(-3.1415926*89.0/180.)


c..  orignal coordinates for central body: trailing edge
	x_in1_te  =  ex_te
	z_in1_te  =  zp


c..  twist angle: twist around leading edge
	ang_twist_w1  =  -2.0
      ang_sf		  =  4.0*3.1415926/180.

c.. design section numbers:
	mmb = 6
	eel = el_in1/real(mmb)

c1
	do i = 1,mmb + 1
        num(i) = number_section(1)
	end do

c... linear interpolate for media airfoil:
	do i =1, number_section(1)
	   ey_up   = sect_up_y  (3,i) - sect_up_y  (4,i)
	   ey_down = sect_down_y(3,i) - sect_down_y(4,i)

	   elup    = ey_up  /real(mmb)
	   eldown  = ey_down/real(mmb)

	   do k =0, mmb
	      up_y  (k+1,i) = sect_up_y  (3,i)  - elup  *real(k)
	      down_y(k+1,i) = sect_down_y(3,i)  - eldown*real(k)
	   end do

	end do

	do k = 0, mmb
	  ez = real(k)*eel 
	  ex_le = x_in1_le	+ ez/(ang_le)
	  ex_te	= x_in1_te	+ ez/(ang_te)
	  eell	= abs(ex_te - ex_le)

c... twist:
	zp = ez + z_in1_le
	em = (zp-0.0)/(el_all-0.0)
        pai = 3.1415926/180.
	ang_twist = (twist_wing(em)*ang_twist_w1)*pai

c	ang_twist = twist_wing(em)*3.1415926/180.*ang_twist_w1
        
c        ang_twist = 0.0    
c        ang_sf=0.0
c... shang fan:
	  d_y = (zp-650.0)*ang_sf

	  do i=1,num(k+1)

	     xm	 = 	rootx(i)*eell
	     ym	 =  up_y(k+1,i)*eell

	     call twist(xm,ym,ang_twist,xe,ye)

	     xx(k+1,i) = xe + ex_le
	     yy(k+1,i) = ye + d_y
	     zz(k+1,i) = ez + z_in1_le
	  end do

	end do

	call face_create_skin(mmb+1,num,xx,yy,zz,nfile)

	do k = 0, mmb
	  ez = real(k)*eel
	  ex_le = x_in1_le + ez/(ang_le)
	  ex_te	= x_in1_te + ez/(ang_te)
	  eell	= abs(ex_te - ex_le)

c... twist:
	zp = ez + z_in1_le
	em = (zp-0.0)/(el_all-0.0)
        pai = 3.1415926/180.
	ang_twist = (twist_wing(em)*ang_twist_w1)*pai

c	ang_twist = twist_wing(em)*3.1415926/180.*ang_twist_w1
c        ang_twist = 0.0    
c        ang_sf=0.0
        
c... shang fan:
	  d_y = (zp-650.0)*ang_sf

	  do i=1,num(k+1)

	     xm	 = 	rootx(i)*eell
	     ym	 =  down_y(k+1,i)*eell

	     call twist(xm,ym,ang_twist,xe,ye)

	     xx(k+1,i) = xe + ex_le
	     yy(k+1,i) = ye	+ d_y
	     zz(k+1,i) = ez + z_in1_le
	  end do

	end do

	call face_create_skin(mmb+1,num,xx,yy,zz,nfile)


111	continue



	close(nfile)
	end

c------------------------------------------------------------------c
      subroutine face_create_skin(n,num,xx,yy,zz,nfile)
c------------------------------------------------------------------c
	dimension xx(400,1000),yy(400,1000),zz(400,1000)
	integer   num(400)

      common /n_vertex/number_vertex
      common /n_edge  /number_edge
      common /n_face  /number_face

	character*25 point
	character*17 edge_nurbs
	character*11 edge_nurbs_end

	character*17 face_head
	character*11 face_d

c---- vetex:
	point         ='vertex create coordinates '

c---- edges:
	edge_nurbs    ='edge create nurbs '
	edge_nurbs_end='interpolate '

c---- faces:
	face_head     ='face create skin '
	face_d	  ='directions '

	npoint0 = number_vertex
	nedge0  = number_edge
	nface0  = number_face

c... edge:
	do m = 1,n
	   nsect = num(m)
	   do i=1, nsect
	      write(nfile,1000)point,xx(m,i),yy(m,i),zz(m,i)
	      number_vertex = number_vertex + 1
	   end do
	   write(nfile,*)edge_nurbs,'\'
	   npoint = number_vertex - nsect
	   do i=npoint+1,npoint+nsect
            call wr_vertex(i,nfile)
	   end do
	   write(nfile,*)edge_nurbs_end
	   number_edge = number_edge + 1
	end do

c... face using skin:
	write(nfile,*)face_head,'\'
	do m = 1,n
	   call wr_edge(nedge0+m,nfile) 
	end do
	write(nfile,*)'directions ','\'
	write(nfile,*)(0,i=1,n)
     
      number_face = number_face + 1
c... generate additional four edges, two edges is none, another is useful. 
	number_edge = number_edge + 4
1000	 format(1x,a25,3f14.5)

	end

c------------------------------------------------------------------c
      subroutine face_create_net(n,num,xx,yy,zz,nfile)
c------------------------------------------------------------------c
	dimension xx(400,1000),yy(400,1000),zz(400,1000)
	integer   num(400)

	dimension x1(1000),y1(1000),z1(1000)
	dimension x2(1000),y2(1000),z2(1000)

      common /n_vertex/number_vertex
      common /n_edge  /number_edge
      common /n_face  /number_face

	character*25 point
	character*17 edge_nurbs
	character*11 edge_nurbs_end

	character*12 face_head
	character*7  face_uedge
	character*7  face_vedge
	character*12 face_u_d
	character*12 face_v_d
	character*4  face_net

c---- vetex:
	point         ='vertex create coordinates '

c---- edges:
	edge_nurbs    ='edge create nurbs '
	edge_nurbs_end='interpolate '

c---- faces:
	face_head     ='face create '
	face_uedge	  ='uedges '
	face_vedge	  ='vedges '
	face_u_d	  ='udirections '
	face_v_d	  ='vdirections	'
	face_net	  ='net '

	npoint0 = number_vertex
	nedge0  = number_edge
	nface0  = number_face

	nsect1 = num(1)
	nsect2 = num(2)

	do i=1,nsect1
	   x1(i) = xx(1,i)
	   y1(i) = yy(1,i)
	   z1(i) = zz(1,i)
	end do

	do i=1,nsect2
	   x2(i) = xx(2,i)
	   y2(i) = yy(2,i)
	   z2(i) = zz(2,i)
	end do

c... edge.1:
	do i=1,nsect1

	write(nfile,1000)point,x1(i),y1(i),z1(i)
	number_vertex = number_vertex + 1

	end do

	write(nfile,*)edge_nurbs,'\'
	do i=npoint0+1,npoint0+nsect1

      call wr_vertex(i,nfile)

	end do
	write(nfile,*)edge_nurbs_end
	number_edge = number_edge + 1

c... edg.2:
	do i=1,nsect2

	write(nfile,1000)point,x2(i),y2(i),z2(i)
	number_vertex = number_vertex + 1

      end do

	write(nfile,*)edge_nurbs,'\'
	do i=npoint0+nsect1+1,npoint0+nsect1+nsect2

      call wr_vertex(i,nfile)

	end do
	write(nfile,*)edge_nurbs_end
	number_edge = number_edge + 1

c...edge.3 j1\j2 is the points for line:
	write(6,*)edge_nurbs,'\'

	j1 = npoint0+1
      call wr_vertex(j1,nfile)

	j2 = npoint0+nsect1 + 1
      call wr_vertex(j2,nfile)

	write(nfile,*)edge_nurbs_end
	number_edge = number_edge + 1

c...edge.4 j3\j4 is the points for line:
	write(6,*)edge_nurbs,'\'

	j3 =npoint0 + nsect1

      call wr_vertex(j3,nfile)

	j4 =npoint0 + nsect1 + nsect2
      call wr_vertex(j4,nfile)

	write(nfile,*)edge_nurbs_end
	number_edge = number_edge + 1

c... face using net:
	write(nfile,*)face_head,face_uedge,'\'
	call wr_edge(nedge0+1,nfile) 
	call wr_edge(nedge0+2,nfile)
	write(nfile,*)'udirections 0 0 ','\'

	write(nfile,*)face_vedge,'\'
	call wr_edge(nedge0+3,nfile) 
	call wr_edge(nedge0+4,nfile)
	write(nfile,*)'vdirections 0 0 ','\'
	write(nfile,*)face_net
	number_face = number_face + 1

1000	 format(1x,a25,3f14.5)

	end




c----------------------------------------------------c
	subroutine wr_vertex(i,nfile)
c----------------------------------------------------c

	if(i.lt.10) then
	write(nfile,100)'"vertex.',i,'"','\'
	end if

	if(i.ge.10.and.i.lt.100) then
	write(nfile,200)'"vertex.',i,'"','\'
	end if

	if(i.ge.100.and.i.lt.1000) then
	write(nfile,300)'"vertex.',i,'"','\'
	end if

	if(i.ge.1000.and.i.lt.10000) then
	write(nfile,400)'"vertex.',i,'"','\'
	end if

	if(i.ge.10000.and.i.lt.100000) then
	write(nfile,500)'"vertex.',i,'"','\'
	end if


100    format(3x,a8,i1,a1,2x,a1)
200    format(3x,a8,i2,a1,2x,a1)
300    format(3x,a8,i3,a1,2x,a1)
400    format(3x,a8,i4,a1,2x,a1)
500    format(3x,a8,i5,a1,2x,a1)
	 end

c----------------------------------------------------c
	subroutine wr_edge(i,nfile)
c----------------------------------------------------c

	if(i.lt.10) then
	write(nfile,100)'"edge.',i,'"','\'
	end if

	if(i.ge.10.and.i.lt.100) then
	write(nfile,200)'"edge.',i,'"','\'
	end if

	if(i.ge.100.and.i.lt.1000) then
	write(nfile,300)'"edge.',i,'"','\'
	end if

	if(i.ge.1000.and.i.lt.10000) then
	write(nfile,400)'"edge.',i,'"','\'
	end if

	if(i.ge.10000.and.i.lt.100000) then
	write(nfile,500)'"edge.',i,'"','\'
	end if


100    format(3x,a6,i1,a1,2x,a1)
200    format(3x,a6,i2,a1,2x,a1)
300    format(3x,a6,i3,a1,2x,a1)
400    format(3x,a6,i4,a1,2x,a1)
500    format(3x,a6,i5,a1,2x,a1)
	 end

c-------------------------------------------------------c
       function twist_body(a)
c-------------------------------------------------------c
	pi2 = 3.1415926*90./180.

      twist_body =0.5-0.5*sin(pi2*(2.*a+1.0)) 
	end

c-------------------------------------------------------c
       function twist_outw(a)
c-------------------------------------------------------c
      twist_outw = 1.+a*.5
	end


c-------------------------------------------------------c
       function twist_wing(a)
c-------------------------------------------------------c
      twist_wing =1.-a
	end



c-------------------------------------------------------c
	 subroutine twist(xm,ym,ang_twist,xe,ye)
c-------------------------------------------------------c
	x0 = 0.0
	y0 = 0.0

	el = sqrt(xm**2. + ym**2.)
	if(el.lt.0.0000000001) then
	   xe = xm
	   ye = ym
	   return
	end if

	cos1 = xm/el
	sin1 = ym/el

	cos2 = cos(ang_twist)
	sin2 = sin(ang_twist)

	cosn = cos1*cos2 - sin1*sin2
	sinn = cos1*sin2 + sin1*cos2

	xe = el*cosn
	ye = el*sinn
c	print*,ym,ye

	end

