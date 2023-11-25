! f2py -c -m river_patch river_patch.F90 --fcompiler=gnu95
! python -m numpy.f2py -c -m river_patch river_patch.F90
!***************************************************************************
subroutine patchms(ix,iy,ngrids,nextX,nextY,uparea,nxtdst,&
    & nx,ny,xlist,ylist,dlist)
    ! get the near by mainstream pixels
    ! n upstream pixels and n number of downstream pixels (2*n+1)
    ! xlist, ylist constist of x, y coordinates
    ! dlist consist of lag distance to x y for ix iy
    ! ix, iy from nearest to farthest 
    !------------------------------------------------
    implicit none
    integer,intent(IN)                           :: ix,iy
    integer,intent(IN)                           :: nx, ny
    integer,intent(IN)                           :: ngrids
    integer,dimension(nx,ny),intent(IN)          :: nextX,nextY !,rivseq
    real,dimension(nx,ny),intent(IN)             :: uparea, nxtdst
    integer,dimension(2*ngrids+1),intent(OUT)    :: xlist, ylist
    real,dimension(2*ngrids+1),intent(OUT)       :: dlist
    !----
    character(len=8)                             :: llon,llat
    character(len=128)                           :: fname
    real,dimension(nx,ny)                        :: weightage
    integer                                      :: info, iix, iiy, fn, k, ios 
    character*128                                :: outdir
    integer,dimension(ngrids)                    :: lx1,ly1,lx2,ly2
    real                                         :: lag_dist
    !===========
    xlist(1)=ix
    ylist(1)=iy
    dlist(1)=0.0
    !===========
    ! upstream girds
    call upstream_pixels(ix,iy,ngrids,nx,ny,nextX,nextY,uparea,lx1,ly1,info)
    ! print*, ix,iy,info
    do k=1, ngrids
        xlist(2*k+1)=lx1(k)
        ylist(2*k+1)=ly1(k)
        call lag_distance(lx1(k),ly1(k),ix,iy,nx,ny,nextX,nextY,nxtdst,lag_dist)
        dlist(2*k+1)=lag_dist
        ! print*, lx1(k), ly1(k), lag_dist
    end do
    !==========
    ! downstream grids
    call downstream_pixels(ix,iy,ngrids,nx,ny,nextX,nextY,lx2,ly2,info)
    ! print*, ix,iy,info
    do k=1, ngrids
        xlist(2*k)=lx2(k)
        ylist(2*k)=ly2(k)
        call lag_distance(lx2(k),ly2(k),ix,iy,nx,ny,nextX,nextY,nxtdst,lag_dist)
        dlist(2*k)=lag_dist
        ! print*, lx2(k), ly2(k), lag_dist
    end do
    ! print*, xlist, ylist, dlist
end subroutine patchms
!*****************************************************************
subroutine upstream(i,j,nx,ny,nextX,nextY,uparea,x,y)
    ! find the upstream pixel with closest upstream area to the i,j
    !--
    implicit none 
    integer,intent(IN)                        :: i,j,nx,ny
    integer,dimension(nx,ny),intent(IN)       :: nextX,nextY !,rivseq
    real,dimension(nx,ny),intent(IN)          :: uparea
    integer,intent(OUT)                       :: x,y
    !--
    real                                      :: dA ! area differnce nextdst,
    integer                                   :: ix,iy,iix,iiy,tx,ty,ud,d
    real                                      :: length,rl
    !--
    x=-9999
    y=-9999
    d=100 ! look at 100*100 box
    dA=1.0e20 ! area differnce
    !--
    !write(*,*)i,j
    do tx=i-d,i+d
        do ty=j-d,j+d
        !write(*,*)tx,ty
        call ixy2iixy(tx,ty, nx, ny, ix, iy)
        !write(*,*)nextX(ix,iy),nextY(ix,iy),ix,iy,uparea(ix,iy),rivseq(ix,iy)
        if (nextX(ix,iy) == i .and. nextY(ix,iy) == j) then
            !write(*,*)ix,iy
            if (abs(uparea(i,j)-uparea(ix,iy)) < dA) then
                dA=abs(uparea(i,j)-uparea(ix,iy))
                x=ix
                y=iy
                ! print*, x,y
            end if
        end if
        end do
    end do
    return
end subroutine upstream
!******************************************************************
subroutine upstream_pixels(i,j,k,nx,ny,nextX,nextY,uparea,lx,ly,info)
    ! list of k number of upstream grids in the mainstream
    implicit none
    !--
    integer                             :: i,j,k,nx,ny
    integer,dimension(nx,ny)            :: nextX,nextY
    real,dimension(nx,ny)               :: uparea
    !--
    !real                                :: threshold,lag_dist
    integer                             :: num,iix,iiy,ix,iy,info
    !real                                :: length,rl
    integer,dimension(k)                :: lx,ly
    !--
    num=1
    iix=i
    iiy=j
    do while (num<=k)
        call upstream(iix,iiy,nx,ny,nextX,nextY,uparea,ix,iy)
        lx(num)=ix
        ly(num)=iy
        iix=ix
        iiy=iy
        num=num+1
        ! print*, iix, iiy, num
    end do
    num=num-1
    info=num
    return
end subroutine upstream_pixels
!******************************************************************
subroutine downstream_pixels(i,j,k,nx,ny,nextX,nextY,lx,ly,info)
    ! k number of downstream grids
    implicit none
    !--
    integer                             :: i,j,k,nx,ny
    integer,dimension(nx,ny)            :: nextX,nextY
    !real,dimension(nx,ny)               :: weight !nextdst
    !--
    !real                                :: threshold,lag_dist
    integer                             :: num,iix,iiy,ix,iy,info
    !eal                                :: length,rl
    integer,dimension(k)                :: lx,ly
    !--
    num=1
    ix=i
    iy=j
    ! lx(num)=ix
    ! ly(num)=iy
    ! num=num+1
    !--
    do while (num<=k)
        iix=ix
        iiy=iy
        ix=nextX(iix,iiy)
        iy=nextY(iix,iiy)
        !write(*,*)ix,iy
        if (ix==-9 .or. iy==-9) then
            num=num+1
            exit
        end if
        if (ix==-10 .or. iy==-10) then
            num=num+1
            exit
        end if
        if (ix==-9999 .or. iy==-9999) then
            num=num+1
            exit
        end if
        lx(num)=ix
        ly(num)=iy
        num=num+1
        ! print*, ix,iy,num
    end do
    num=num-1
    info=num
    return
end subroutine downstream_pixels
!*****************************************************************
function roundx(ix, nx)
    implicit none
    !-- for input -----------
    integer                     ix, nx
    !-- for output ----------
    integer                     roundx
    !------------------------
    if (ix .ge. 1) then
      roundx = ix - int((ix -1)/nx)*nx
    else
      roundx = nx - abs(mod(ix,nx))
    end if 
    return
end function roundx
!*****************************************************************
subroutine ixy2iixy(ix,iy, nx, ny, iix, iiy)
    implicit none
    !- for input -----------------
    integer                   ix, iy, nx, ny
    !- for output ----------------
    integer                   iix, iiy,roundx
    !-----------------------------
    if (iy .lt. 1) then
      iiy = 2 - iy
      iix = ix + int(nx/2.0)
      iix = roundx(iix, nx)
    else if (iy .gt. ny) then
      iiy = 2*ny -iy
      iix = ix + int(nx/2.0)
      iix = roundx(iix, nx)
    else
      iiy = iy
      iix = roundx(ix, nx)
    end if
    return
end subroutine ixy2iixy
!*****************************************************************
subroutine lag_distance(i,j,x,y,nx,ny,nextX,nextY,nextdst,lag_dist)
    implicit none 
    !--
    integer                             :: i,j,x,y,nx,ny
    integer,dimension(nx,ny)            :: nextX,nextY
    real,dimension(nx,ny)               :: nextdst
    !--
    real                                :: lag_dist
    integer                             :: ix,iy,iix,iiy,tx,ty,ud
    real                                :: length,rl
    !--
    if (i==x .and. j==y) then
      ud=0
    else
      ud=-1
    end if
    !--
    if (ud==-1) then
      tx=x
      ty=y
      ix=i
      iy=j
      length=0.0
      lag_dist=0.0
      do while (ix/=tx .or. iy/=ty) 
        iix=ix
        iiy=iy 
        ix=nextX(iix,iiy)
        iy=nextY(iix,iiy)
        if (ix==-9 .or. iy==-9) then ! river mouth to the sea
          ud=+1
          exit
        end if
        if (ix==-10 .or. iy==-10) then ! inland termination
          ud=+1
          exit
        end if
        if (ix==-9999 .or. iy==-9999) then
          ud=+1
          exit
        end if
        !-- distance to the next grid
        !print*, nextdst(ix,iy)
        rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
        length=length+rl!/2.0
      end do
    end if
    !--
    if (ud==+1) then
      tx=i
      ty=j
      ix=x
      iy=y
      length=0.0
      do while (ix/=tx .or. iy/=ty) 
        iix=ix
        iiy=iy
        ix=nextX(iix,iiy)
        iy=nextY(iix,iiy)
        !---
        if (ix==-9 .or. iy==-9) then ! river mouth to the sea
          ud=-9999
          exit
        end if
        if (ix==-10 .or. iy==-10) then ! inland termination
          ud=+1
          exit
        end if
        if (ix==-9999 .or. iy==-9999) then
          ud=-9999
          exit
        end if
        !-- half of the present grid
        rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
        length=length+rl!/2.0
      end do
    end if
    !-- 
    if (ud==-9999) then
      lag_dist=-9999
    elseif (ud==0) then
      lag_dist=0.0
    else
      lag_dist=length
    end if
    !---
    return
    !---
end subroutine lag_distance
!*****************************************************************
