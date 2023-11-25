program read_dis
    !==========================================
    implicit none
    ! index CaMa
    integer                       ::  iXX, iYY !, jXX, jYY, kXX, kYY, uXX, uYY
    integer                       ::  nXX, nYY, nFL                   !! x-y matrix GLOBAL
    real                          ::  gsize, csize                    !! grid size [degree]
    real                          ::  west, north, east, south
    character*128                 ::  buf
    character*128                 ::  camadir, map, tag               !! CaMa dir , map[e.g.:glb_15min], tag[e.g. 1min, 15sec, 3sec]
    character*128                 ::  river, station, id
    integer                       ::  ix1, iy1, ix2, iy2
    ! integer                       ::  ngrid
    integer                       ::  ios
    integer                       ::  cnum, num, mwin

    integer                       ::  year, mon, day
    integer                       ::  syear, smon, sday, eyear, emon, eday
    character*128                 ::  indir, outdir, expname, mapname, fname

    integer                       ::  isleap, leap, monthday, ens_num, ens
    integer                       ::  N, i
    character*8                   ::  yyyymmdd
    character*8,allocatable       ::  daylist(:)
    character*3                   ::  numch

    real,allocatable              ::  opn(:,:,:), opnout(:,:,:), asm(:,:,:), asmout(:,:,:)
    character*128                 ::  rlist, fout
    character*128                 ::  buffer


    call getarg(1,buf)
    read(buf,"(A)") expname
    call getarg(2,buf)
    read(buf,"(A)") mapname
    call getarg(3,buf)
    read(buf,*) syear
    call getarg(4,buf)
    read(buf,*) smon
    call getarg(5,buf)
    read(buf,*) sday
    call getarg(6,buf)
    read(buf,*) eyear
    call getarg(7,buf)
    read(buf,*) emon
    call getarg(8,buf)
    read(buf,*) eday
    call getarg(9,buf)
    read(buf,*) ens_num ! number of ensemble
    call getarg(10,buf)
    read(buf,*) N
    call getarg(11,buf)
    read(buf,"(A)") camadir
    call getarg(12,buf)
    read(buf,"(A)") indir
    call getarg(13,buf)
    read(buf,"(A)") outdir
    !==
    fname=trim(camadir)//"/map/"//trim(mapname)//"/params.txt"
    print *, fname
    open(11,file=fname,form='formatted')
    read(11,*) nXX
    read(11,*) nYY
    read(11,*) nFL
    read(11,*) gsize
    read(11,*) west
    read(11,*) east
    read(11,*) south
    read(11,*) north
    close(11)
    !==
    allocate(opn(nXX,nYY,ens_num),opnout(N,nXX,nYY),asm(nXX,nYY,ens_num),asmout(N,nXX,nYY),daylist(N))
    !================
    num=1
    do year=syear,eyear
        isleap=leap(year)
        do mon=smon, emon
            eday=monthday(mon,isleap)
            do day=1, eday
                write(yyyymmdd,'(i4.4,i2.2,i2.2)') year,mon,day
                do ens=1, ens_num
                    write(numch,'(i3.3)') ens
                    ! courrpted
                    fname=trim(indir)//"/"//trim(expname)//"/assim_out/outflw/open/outflw"//trim(yyyymmdd)//"_"//numch//".bin"
                    open(34,file=fname,form="unformatted",access="direct",recl=4*nXX*nYY,status="old",iostat=ios)
                    if(ios==0)then
                        read(34,rec=1) opn(:,:,ens)
                    else
                        write(*,*) "no: ",trim(fname)
                    end if
                    close(34)
                    ! assimilated
                    fname=trim(indir)//"/"//trim(expname)//"/assim_out/outflw/assim/outflw"//trim(yyyymmdd)//"_"//numch//".bin"
                    open(34,file=fname,form="unformatted",access="direct",recl=4*nXX*nYY,status="old",iostat=ios)
                    if(ios==0)then
                        read(34,rec=1) asm(:,:,ens)
                    else
                        write(*,*) "no: ",trim(fname)
                    end if
                    close(34)
                end do
                opnout(num,:,:)=sum(opn,DIM=3)/(real(ens_num)+1e-20)
                asmout(num,:,:)=sum(asm,DIM=3)/(real(ens_num)+1e-20)
                daylist(num)=trim(yyyymmdd)
                num=num+1
            end do
        end do
    end do
    ! ===============================================
    ! read data 
    ! ===============================================
    rlist=trim(camadir)//"/map/"//trim(mapname)//"/grdc_loc.txt"
    open(11, file=rlist, form='formatted')
    read(11,*)
    !----
1000 continue
    read(11,'(A)',end=1090) buffer
    id=buffer(1:7)
    river=buffer(9:43)
    station=buffer(44:86)
    read(buffer( 87: 94),'(i7)') ix1
    read(buffer( 95:102),'(i7)') iy1
    read(buffer(103:110),'(i7)') ix2
    read(buffer(111:118),'(i7)') iy2
    print*, trim(adjustl(station)), ix1, iy1
    fout=trim(adjustl(outdir))//"/"//trim(expname)//"/outflow/"//trim(adjustl(id))//".txt"
    open(72,file=fout,status='replace')
    do i=1, N 
        if (ix2==-9999 .or. iy2==-9999) then
            ! print'(a8,4x,f10.2,4x,f10.2)', daylist(i), asmout(i,ix1,iy1), opnout(i,ix1,iy1)
            write(72,'(a8,4x,f15.2,4x,f15.2)')daylist(i), asmout(i,ix1,iy1), opnout(i,ix1,iy1)
        else
            ! print'(a8,4x,f10.2,4x,f10.2)', daylist(i), asmout(i,ix1,iy1)+asmout(i,ix2,iy2), opnout(i,ix1,iy1)+opnout(i,ix2,iy2)
            write(72,'(a8,4x,f15.2,4x,f15.2)')daylist(i), asmout(i,ix1,iy1)+asmout(i,ix2,iy2), opnout(i,ix1,iy1)+opnout(i,ix2,iy2)
        end if
    end do
    close(72)
    goto 1000
1090 continue
    close(11)
    deallocate(opn,opnout,asm,asmout,daylist)
end program read_dis
!*********************************
function leap(year)
    implicit none
    integer                        :: year,leap
    !real                           :: mod
    !--
    ! days of the year
    ! normal : 365
    ! leap   : 366
    leap=0
    if (mod(dble(year),4.0)   == 0) leap=1
    if (mod(dble(year),100.0) == 0) leap=0
    if (mod(dble(year),400.0) == 0) leap=1
    return
end function leap
!*********************************
function monthday(imon,leap)
    implicit none
    integer                        :: imon, leap, monthday
    ! ================================================
    ! Calculation for months except February
    ! ================================================
    if(imon.eq.4.or.imon.eq.6.or.imon.eq.9.or.imon.eq.11)then
        monthday=30
    else
        monthday=31
    end if
    ! ================================================
    ! Calculation for February
    ! ================================================
    if ( leap.eq.1)then
        if (imon.eq.2)then
            monthday=29
        end if
    else
        if (imon.eq.2)then
            monthday=28
        end if
    end if
    return
end function monthday
!*********************************