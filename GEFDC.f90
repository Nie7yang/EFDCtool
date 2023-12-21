    !****************************************************************************
    !
    !  PROGRAM: GEFDC
    !
    !  PURPOSE:  Updated GEFDC.
    !            1.Adopted free format;
    !            2.Replaced some outdated common syntax;
    !            3.Corrected the boundary issue of CIJ(I, J) in the
    !              difference format processing.
    !
    !****************************************************************************
    !  Updated by Nie Qiyang
    !             Hokkaido University
    !             2023.12.14
    !****************************************************************************


    !****************************************************************************
    !  Module 模块
    !  
    !****************************************************************************

    Module OneModule
    implicit none

    real(kind=8),allocatable :: depd(:), depe(:), depn(:)
    real(kind=8) :: cdep, radm
    end Module OneModule

    Module TwoModule
    implicit none

    integer :: ndepdat
    end Module TwoModule

    Module VegModule
    implicit none

    integer :: nvegdat, nvegtyp

    real(kind=8) ,allocatable:: nvegd(:), vege(:), vegn(:)
    integer ,allocatable:: nnveg(:)
    end Module VegModule

    Module CutmModule
    implicit none

    real(kind=8) :: rk, phi0, cmerid, ae, ecc2, ecc2p
    end Module CutmModule

    Module CsmdModule
    implicit none

    real(kind=8) :: smd
    integer :: ib, ie, jb, je, ijsmd, ismd, jsmd
    end Module CsmdModule

    Module CxyModule
    implicit none

    real(kind=8),allocatable :: x(:,:), y(:,:)

    end Module CxyModule

    Module FuncModule

    CONTAINS

    real(kind=8) function xutmbay(rlong, rlat)
    use CutmModule
    real(kind=8):: rlong, rlat, x ,y, xutm, rn, t, c, a

    x = 500000.  ! 设置原点的东坐标为 500,000 米
    y = 0.       ! 设置原点的北坐标为 0 米
    ! 检查纬度是否在有效范围之外
    if (rlat>80. .Or. rlat<-80.) Then
        xutm = 0. ! 将 UTM 东坐标设置为 0 米
    else
        rn = ae/sqrt(1-ecc2*(sin(rlat*dtr))**2)  ! 计算子午圈的曲率半径
        t = tan(rlat*dtr)**2    ! 计算纬度的正切值的平方
        c = ecc2p*cos(rlat*dtr)**2   ! 计算偏心率的平方乘以纬度余弦的平方
        a = cos(rlat*dtr)*(rlong-cmerid)*dtr  ! 计算从中央子午线的经度差
        xutm = rk*rn*(a+(1-t+c)*(a**3)/6.+(5.-18.*t+(t**2)+72*c-58.*ecc2p)*(a**5)/120.)  ! 计算 UTM 东坐标
    endif
    xutmbay = x - xutm  ! 计算 UTM 东坐标
    return
    end function xutmbay

    real(kind=8) function yutmbay(rlong, rlat)
    use CutmModule
    real(kind=8) rlong, rlat, x ,y, rn, t, c, a, rm, rm0, tempa, tempb, tempc, tempd, tempe, temp1, temp2, temp3, temp

    x = 500000.
    y = 0.
    rn = ae/sqrt(1-ecc2*sin(rlat*dtr)**2)
    t = tan(rlat*dtr)**2
    c = ecc2p*cos(rlat*dtr)**2
    a = cos(rlat*dtr)*(rlong-cmerid)*dtr
    rm = ae*((1-ecc2/4-3*ecc2**2/64-5*ecc2**3/256)*rlat*dtr-(3*ecc2/8+3*ecc2**2/32+45*ecc2**3/1024)*sin(2*rlat*dtr)+(15*ecc2**2/256+45*ecc2**3/1024)*sin(4*rlat*dtr)-(35*ecc2**3/3072)*sin(6*rlat*dtr))
    rm0 = ae*((1-ecc2/4-3*ecc2**2/64-5*ecc2**3/256)*phi0*dtr-(3*ecc2/8+3*ecc2**2/32+45*ecc2**3/1024)*sin(2*phi0*dtr)+(15*ecc2**2/256+45*ecc2**3/1024)*sin(4*phi0*dtr)-(35*ecc2**3/3072)*sin(6*phi0*dtr))
     ! 检查纬度是否在极点
    if (rlat==90. .Or. rlat==-90.) Then
        yutmbay = rk*(rm-rm0)  ! 计算 UTM 北坐标
    else
        tempa = a**6/720.
        tempb = 330.*ecc2p
        tempc = 600.*c
        tempd = t**2
        tempe = 58.*t
        temp1 = (61.-tempe+tempd+tempc-tempb)*tempa
        tempa = a**4/24.
        tempb = 4*c**2
        tempc = 9.*c
        temp2 = (5.-t+tempc+tempb)*tempa
        temp3 = a**2/2
        temp = temp1 + temp2 + temp3
        temp = tan(rlat*dtr)*temp
        temp = rn*temp
        temp = rm + temp
        yutmbay = rk*temp
    endif
    return
    end function yutmbay

    real(kind=8) function fib(yy, j)
    implicit none
    real(kind=8):: yy
    integer:: j
    if (yy>=20.6 .And. yy<=35.6) Then
        fib = 9.*(yy-21.)/15. + 7.
        return
    endif
    write (6,'(" function FIE OUT OF BOUNDS YY,J = ", F10.4, I8/)') yy, j
    return
    end function fib

    real(kind=8) function fie(yy, j)
    implicit none
    real(kind=8):: yy
    integer:: j
    if (yy>=9.1 .And. yy<=23.6) Then
        fie = 31.
        return
    endif
    write (6,'(" function FIE OUT OF BOUNDS YY,J = ", F10.4, I8/)') yy, j
    return
    end function fie

    real(kind=8) function gjb(xx, i)
    implicit none
    real(kind=8):: xx, x, ctmp, dtmp
    integer:: i
    if (xx>=6.76 .And. xx<8.76) Then
        x = xx - 6.76
        gjb = 20.6 - 0.6*x - (2.254-0.203*x)*x*x/7.7
        return
    endif
    if (xx>=8.76 .And. xx<14.7) Then
        gjb = -11.2*(xx-7.)/7.7 + 21.
        return
    endif
    if (xx>=14.7 .And. xx<19.4) Then
        x = xx - 14.7
        ctmp = 6.764968/(4.7*4.7)
        dtmp = -2.028605/(4.7*4.7*4.7)
        gjb = 9.8 - 11.2*x/7.7 + (ctmp+dtmp*x)*x*x
        return
    endif
    if (xx>=19.4 .And. xx<=29.0) Then
        gjb = 1.5*(xx-19.4)/11.6 + 7.7
        return
    endif
    if (xx>=29.0 .And. xx<=31.) Then
        x = xx - 31.
        gjb = 9.1 - (0.63+0.085*x)*x*x/11.6
        return
    endif
    write (6, '(" function GJE OUT OF BOUNDS XX,I = ", F10.4, I8/)') xx, i
    return
    end function gjb

    real(kind=8) function gje(xx, i)
    implicit none
    real(kind=8):: xx, x
    integer:: i
    if (xx>=15.76 .And. xx<17.76) Then
        x = xx - 15.76
        gje = 35.6 - 0.6*x - (1.696-0.082*x)*x*x/7.5
        return
    endif
    if (xx>=17.76 .And. xx<22.5) Then
        gje = -10.3*(xx-16.)/7.5 + 36.
        return
    endif
    if (xx>=22.5 .And. xx<24.5) Then
        x = xx - 22.5
        gje = (203.05-10.3*x+2.*x*x)/7.5
        return
    endif
    if (xx>=24.5 .And. xx<29.0) Then
        gje = -2.3*(xx-23.5)/7.5 + 25.7
        return
    endif
    if (xx>=29.0 .And. xx<=31.0) Then
        x = xx - 31.
        gje = 23.6 + (1.175+0.2*x)*x*x/7.5
        return
    endif
    write (6, '(" function GJE OUT OF BOUNDS XX,I = ", F10.4, I8/ )') xx, i
    return
    end function gje

    end Module FuncModule

    !****************************************************************************
    !  �����򲿷�
    !
    !****************************************************************************
    Program gefdc

    use FuncModule
    use OneModule
    use TwoModule
    use VegModule
    use CutmModule
    use CsmdModule
    use CxyModule

    parameter ( imind=1, jmind=1, imaxd=150, jmaxd=150, iggdim=150, jggdim=150, nbpd=1000, ngim=2500, nwcdim=2500)
    integer(kind=8), parameter :: nddmax=225000, nvdmax=1000
    parameter ( pii=3.1415926535898, dtr=0.017453292519943)
    real(kind=8), dimension(imind:imaxd, jmind:jmaxd) :: xn, yn, hi, hj, rki, rkj, rkii, rkjj, rsgi, cij, nvegij
    integer, dimension(ngim) :: ired, jred, iblk, jblk
    integer, dimension(imind:imaxd, jmind:jmaxd) :: ksgi, ksbp, ijct
    real(kind=8), dimension(imind:imaxd, jmind:jmaxd) :: xcell, dlondd, ycell, dlatdd, depfix, depcc
    integer, dimension(nbpd) :: ibp, jbp
    integer, dimension(nwcdim) :: icomp, jcomp
    integer, dimension(iggdim, jggdim) :: ijctg
    real(kind=8), dimension(nwcdim) :: xcomp, ycomp

    integer:: i,j,k,l,m,n, ntype, nbpp, imin, imax, jmin, jmax, ic, jc, isgg, igm, jgm, nwtgg, itrxm, itrhm, itrkm, itrgm, ndepsm, ndepsmf, isirki, jsirki, isihihj, jsihihj, n7relax, nxyit, itn7max, isidep, isidptyp, isveg, ilt, jlt, ift, jft
    real(kind=8):: dxcg, dycg, cdlon1, cdlon2, cdlon3, cdlat1, cdlat2, cdlat3, depmin, ddatadj, rpx, rpk, rph, rsqxm, rsqkm, rsqkim, rsqhm, rsqhim, rsqhjm, xshift, yshift, hscale, rkjdki, angoro, rp7, serrmax,xibjb, yibjb,   surfelv, ytmp, xtmp, xlnutme, yltutmn, radsq1, radsq2, nvtmp
    character(len=80):: title, dummy
    integer:: handle_read_error

    allocate ( depd(nddmax), depe(nddmax), depn(nddmax), nvegd(nvdmax), vege(nvdmax), vegn(nvdmax),nnveg(0:12) )
    allocate ( x(imind:imaxd, jmind:jmaxd), y(imind:imaxd, jmind:jmaxd) )

    open (2, file='depdat.inp', status='unknown')
    open (3, file='cell.inp', status='unknown')
    open (4, file='gefdc.inp', status='unknown')
    open (7, file='gefdc.out', status='unknown')
    open (8, file='grid.ixy', status='unknown')
    open (9, file='grid.jxy', status='unknown')
    open (10, file='grid.mask', status='unknown')
    open (12, file='grid.cord', status='unknown')
    open (13, file='grid.init', status='unknown')
    open (14, file='dxdy.out', status='unknown')
    open (15, file='dxdy.diag', status='unknown')
    open (16, file='lxly.out', status='unknown')
    open (66, file='gefdc.log', status='unknown')
    open (67, file='depint.log', status='unknown')
    open (25, file='grid.dxf', status='unknown')
    open (26, file='init.dxf', status='unknown')
    open (27, file='data.plt', status='unknown')
    open (89, file='depspc.out', status='unknown')
    write (25, "('  0', /, 'SECTION', /, '  2', /, 'ENTITIES')")
    write (26, "('  0', /, 'SECTION', /, '  2', /, 'ENTITIES')")
    read (4, *) dummy
    read (4, *) dummy
    read (4, *) title
    read (4, *) dummy
    read (4, *) dummy
    read (4, *) ntype, nbpp, imin, imax, jmin, jmax, ic, jc
    write (7, *) ntype, nbpp, imin, imax, jmin, jmax, ic, jc
    read (4, *) dummy
    read (4, *) dummy
    read (4, *) isgg, igm, jgm, dxcg, dycg, nwtgg
    write (7, *) isgg, igm, jgm, dxcg, dycg, nwtgg
    read (4, *) dummy
    read (4, *) dummy
    read (4, *) cdlon1, cdlon2, cdlon3, cdlat1, cdlat2, cdlat3
    write (7, *) cdlon1, cdlon2, cdlon3, cdlat1, cdlat2, cdlat3
    read (4, *) dummy
    read (4, *) dummy
    read (4, *) itrxm, itrhm, itrkm, itrgm, ndepsm, ndepsmf, depmin, ddatadj
    write (7, *) itrxm, itrhm, itrkm, itrgm, ndepsm, ndepsmf, depmin, ddatadj
    read (4, *) dummy
    read (4, *) dummy
    read (4, *) rpx, rpk, rph, rsqxm, rsqkm, rsqkim, rsqhm, rsqhim, rsqhjm
    write (7, *) rpx, rpk, rph, rsqxm, rsqkm, rsqkim, rsqhm, rsqhim, rsqhjm
    read (4, *) dummy
    read (4, *) dummy
    read (4, *) xshift, yshift, hscale, rkjdki, angoro
    write (7, *) xshift, yshift, hscale, rkjdki, angoro
    read (4, *) dummy
    read (4, *) dummy
    read (4, *) isirki, jsirki, isihihj, jsihihj
    write (7, *) isirki, jsirki, isihihj, jsihihj
    read (4, *) dummy
    read (4, *) dummy
    if (ntype==7) Then
        read (4, *) ib, ie, jb, je, n7relax, nxyit, itn7max, ijsmd, ismd, jsmd, rp7, serrmax
    endif
    if (ijsmd/=0) Then
        ismd = 0
        jsmd = 0
    endif
    if (ismd/=0) Then
        ijsmd = 0
        jsmd = 0
    endif
    if (jsmd/=0) Then
        ijsmd = 0
        ismd = 0
    endif
    read (4, *) dummy
    read (4, *) dummy
    if (ntype==7) Then
        read (4, *) xibjb, yibjb
        read (4, *) xiejb, yiejb
        read (4, *) xieje, yieje
        read (4, *) xibje, yibje
    endif
    read (4, *) dummy
    read (4, *) dummy
    read (4, *) isidep, ndepdat, cdep, radm, isidptyp, surfelv, isveg, nvegdat, nvegtyp
    write (7, *) isidep, ndepdat, cdep, radm, isidptyp, surfelv, isveg, nvegdat, nvegtyp
    imaxo = imax
    jmaxo = jmax
    imino = imin
    jmino = jmin
    if (isidep==1) Then
        ncount = 0
        do n = 1, ndepdat
            ncount = ncount + 1
            read (2, *, iostat=handle_read_error) depe(n), depn(n), deptmp
            if (handle_read_error/= 0)then
                write(6, "('read ERROR ON file depdat.inp AT NCOUNT =', I10//)")
                stop
            endif
            depd(n) = ddatadj + deptmp
            if (isidptyp==1) depd(n) = max(depmin, depd(n))
            if (isidptyp==3) depd(n) = -depd(n)
        enddo
        write (6, "(I10, ' DEPTH DATA POINTS read IN', //)") ncount
    endif

    if (isveg==1) Then
        open (90, file='vegdat.inp', status='unknown')
        ncount = 0
        do n = 1, nvegdat
            ncount = ncount + 1
            read (90, *, iostat=handle_read_error) vege(n), vegn(n), nvegd(n)
            if (handle_read_error/= 0)then
                write(6, "('read ERROR ON file vegdat.inp AT NCOUNT =', I10//)")
                stop
            endif
        enddo
        close (90)
    endif

    if (isveg==1) write (6, "(I10, ' VEGATATION DATA POINTS read IN', //)") ncount
    if (ntype>=5 .And. ntype<=6) Then
        nbp = nbpp
    else
        nbp = 2*nbpp - 2
    endif
    if (ntype==1) Then
        imax = 2*imax - imin
    endif
    if (ntype==2) Then
        jmax = 2*jmax - jmin
    endif
    if (ntype==3) Then
        imin = 2*imin - imax
    endif
    if (ntype==4) Then
        jmin = 2*jmin - jmax
    endif
    do j = jmin, jmax
        do i = imin, imax
            x(i, j) = 0.
            y(i, j) = 0.
            xn(i, j) = 0.
            yn(i, j) = 0.
            hi(i, j) = 1.
            hj(i, j) = 1.
            rki(i, j) = 1.
            rkii(i, j) = 1.
            rkj(i, j) = rkjdki
            rkjj(i, j) = rkjdki
            ksbp(i, j) = 0
            ksgi(i, j) = 0
            rsgi(i, j) = 0.
            depcc(i, j) = 0.
            depfix(i, j) = 0.
            nvegij(i, j) = 0
        enddo
    enddo

    !********** (ntype>=8) if **********
    if (ntype>=8) then
        do is = 1, 4
            read (3, "(130X)")
        enddo
        read (3, "(I3, 2X, 125I1)") jctmp
        close (3)
        open (3, file='cell.inp', status='unknown')
        do is = 1, 4
            read (3, "(130X)")
        enddo
        if (jctmp/=jc) Then
            do jt = 1, jc, 120
                jf = jt
                jlast = jt + 119
                if (jlast>jc) jlast = jc
                write (7, "(1X, '  CELL TYPE ARRAY,J=', I5, 2X, 'TO J=', I5, //)") jf, jlast
                do i = 1, ic
                    read (3, "(120I1)")(ijct(i,j), j=jf, jlast)
                    write (7, "(1X, I3, 2X, 125I1)") i, (ijct(i,j), j=jf, jlast)
                enddo
            enddo
        else
            do it = 1, ic, 125
                ifirst = it
                ilast = it + 124
                if (ilast>ic) ilast = ic
                write (7, "(1X, '  CELL TYPE ARRAY,I=', I5, 2X, 'TO I=', I5, //)") ifirst, ilast
                do j = jc, 1, -1
                    read (3, "(I3, 2X, 125I1)") jdumy, (ijct(i,j), i=ifirst, ilast)
                    write (7, "(1X, I3, 2X, 125I1)") jdumy, (ijct(i,j), i=ifirst, ilast)
                enddo
            enddo
        endif
        open (18, file='gcellmap.inp', status='unknown')
        open (30, file='salt.inp', status='unknown')
        nwtgg = 1
        write (14, "('C dxdy.inp file, in free format across columns')")
        write (14, "('C')")
        write (14, "('C     I     J        DX            DY           ', 1X, 'DEPTH     BOTTOM ELEV      ZROUGH  VEG TYPE')")
        write (14, "('C')")
        write (16, "('C lxly.inp file, in free format across line')")
        write (16, "('C')")
        write (16, "('C    I     J    XLNUTME       YLTUTMN        CCUE', 1X, '           CCVE          CCUN         CCVN')")
        write (16, "('C')")
        write (18, "('C gcellmap.inp file, in free format across columns')")
        write (18, "('C')")
        write (18, "('  IGRAPHIC  JGRAPHIC     ICOMP     JCOMP         RMIN')")
        write (18, "('C')")
        write (18, "(4I10, 5X, F12.4)") igm, jgm, nwtgg
        if (ntype/=9) Then
            n999 = 0
            depmax = 0.
            nwcells = 0.
            do j = 1, jc
                do i = 1, ic
                    if (ijct(i,j)>=1 .And. ijct(i,j)<=5) Then
                        nwcells = nwcells + 1
                        ceu = 1.0
                        cev = 0.0
                        cnu = 0.0
                        cnv = 1.0
                        xlnutme = cdlon1 + (cdlon2*float(i)+cdlon3)/60.
                        yltutmn = cdlat1 + (cdlat2*float(j)+cdlat3)/60.
                        radsq1 = dxcg*dycg
                        radsq2 = dxcg*dycg
                        if (isidep==1) Then
                            Call raddep(i, j, xlnutme, yltutmn, radsq1, radsq2, depth)
                            if (depth==-999.) Then
                                n999 = n999 + 1
                                write (67, "('NO DEP DATA AT I,J,X,Y = ', 2I5, 2E14.5)") i, j, xlnutme, ylnutmn
                            endif
                            if (isidptyp==1) Then
                                belv = -1.*depth
                                zrough = 0.0
                                nvegdum = 0
                                depmax = max(depmax, depth)
                            else
                                belv = depth
                                depth = surfelv - belv
                                zrough = 0.0
                                nvegdum = 0
                                depmax = max(depmax, depth)
                            endif
                        else
                            depth = 0.
                            belv = 0.0
                            zrough = 0.0
                            nvegdum = 0
                        endif
                        wndshe = 1.0

                        write (14, "(1X, I5, 2X, I5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, I3)") i, j, dxcg, dycg, depth, belv, zrough, nvegdum
                        write (6, "(1X, I5, 2X, I5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, I3)") i, j, dxcg, dycg, depth, belv, zrough, nvegdum
                        write (66, "(1X, I5, 2X, I5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, I3)") i, j, dxcg, dycg, depth, belv, zrough, nvegdum
                        write (16, "(1X, I5, 1X, I5, 7(1X,E13.6))") i, j, xlnutme, yltutmn, ceu, cev, cnu, cnv, wndshe
                        if (ijct(i,j)==1) Then
                            write (27, "(2F12.3, '  -1')") x(i, j), y(i, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j), y(i+1, j)
                            write (27, "(2F12.3, '   1')") x(i, j+1), y(i, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j), y(i, j)
                        endif
                        if (ijct(i,j)==4) Then
                            write (27, "(2F12.3, '  -1')") x(i, j), y(i, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j), y(i+1, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j+1), y(i+1, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j), y(i, j)
                        endif
                        if (ijct(i,j)==3) Then
                            write (27, "(2F12.3, '  -1')") x(i+1, j), y(i+1, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j+1), y(i+1, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j+1), y(i, j+1)
                            write (27, "(2F12.3, '   1')") x(i+1, j), y(i+1, j)
                        endif
                        if (ijct(i,j)==2) Then
                            write (27, "(2F12.3, '  -1')") x(i, j), y(i, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j+1), y(i+1, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j+1), y(i, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j), y(i, j)
                        endif
                        if (ijct(i,j)==5) Then
                            write (27, "(2F12.3, '  -1')") x(i, j), y(i, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j), y(i+1, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j+1), y(i+1, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j+1), y(i, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j), y(i, j)
                        endif
                        if (ijct(i,j)==1) Then
                            write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==4) Then
                            write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==3) Then
                            write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (25, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==2) Then
                            write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==5) Then
                            write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'SEQEND')")
                        endif
                        icompt = i
                        jcompt = j
                        rmin = 0.
                        write (18, "(4I10, 5X, F12.4)") i, j, icompt, jcompt, rmin
                    endif
                enddo
            enddo
        else
            Call vutmbay
            do j = 1, jc
                do i = 1, ic
                    xlnutme = abs(cdlon1+(cdlon2*float(i)+cdlon3)/60.)
                    yltutmn = cdlat1 + (cdlat2*float(j)+cdlat3)/60.
                    xcell(i, j) = xutmbay(xlnutme, yltutmn)
                    ycell(i, j) = yutmbay(xlnutme, yltutmn)
                    xcell(i, j) = xcell(i, j)*1.E-3
                    ycell(i, j) = (ycell(i,j)-0.4E+7)*1.E-3
                enddo
            enddo
            do j = 2, jc
                do i = 2, ic
                    if (ijct(i,j)==0) Then
                        x(i, j) = 0.
                        y(i, j) = 0.
                    else
                        x(i, j) = 0.25*(xcell(i,j)+xcell(i-1,j)+xcell(i,j-1)+xcell(i-1,j-1))
                        y(i, j) = 0.25*(ycell(i,j)+ycell(i-1,j)+ycell(i,j-1)+ycell(i-1,j-1))
                    endif
                enddo
            enddo
            n999 = 0
            depmax = 0.
            nwcells = 0
            do j = 1, jc
                do i = 1, ic
                    xlnutme = abs(cdlon1+(cdlon2*float(i)+cdlon3)/60.)
                    yltutmn = cdlat1 + (cdlat2*float(j)+cdlat3)/60.
                    dlondd(i, j) = xlnutme
                    dlatdd(i, j) = yltutmn
                    xcell(i, j) = xutmbay(xlnutme, yltutmn)
                    ycell(i, j) = yutmbay(xlnutme, yltutmn)
                    if (ijct(i,j)>=1 .And. ijct(i,j)<=5) nwcells = nwcells + 1
                enddo
            enddo
            do j = 1, jc
                do i = 1, ic
                    if (ijct(i,j)>=1 .And. ijct(i,j)<=5) Then
                        ceu = 1.0
                        cev = 0.0
                        cnu = 0.0
                        cnv = 1.0
                        dxcctr = 0.5*abs(xcell(i+1,j)-xcell(i-1,j))
                        dycctr = 0.5*abs(ycell(i,j+1)-ycell(i,j-1))
                        xlnutme = xcell(i, j)*1.E-3
                        yltutmn = (ycell(i,j)-0.4E+7)*1.E-3
                        radsq1 = dxcctr*dycctr*1.E-6
                        radsq2 = dxcctr*dycctr*1.E-6
                        if (isidep==1) Then
                            Call raddep(i, j, xlnutme, yltutmn, radsq1, radsq2, depth)
                            if (depth==-999.) Then
                                write (67, "(2I5, 2X, F10.4, 2X, F10.4, 2X, F12.4, 2X, F12.4)") i, j, dlondd(i, j), dlatdd(i, j), xlnutme, yltutmn
                            endif
                            if (depth<1.0) depth = 1.0
                            belv = -1.*depth
                            zrough = 0.0
                            nvegdum = 0
                            depmax = max(depmax, depth)
                        else
                            depth = 0.
                            belv = 0.0
                            zrough = 0.0
                            nvegdum = 0
                        endif
                        wndshe = 1.0
                        write (14, "(1X, I5, 2X, I5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, I3)") i, j, dxcctr, dycctr, depth, belv, zrough, nvegdum
                        write (6, "(1X, I5, 2X, I5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, I3)") i, j, dxcctr, dycctr, depth, belv, zrough, nvegdum
                        write (66, "(1X, I5, 2X, I5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, I3)") i, j, dxcctr, dycctr, depth, belv, zrough, nvegdum
                        write (16, "(1X, I5, 1X, I5, 7(1X,E13.6))") i, j, xlnutme, yltutmn, ceu, cev, cnu, cnv, wndshe
                        if (ijct(i,j)==1) Then
                            write (27, "(2F12.3, '  -1')") x(i, j), y(i, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j), y(i+1, j)
                            write (27, "(2F12.3, '   1')") x(i, j+1), y(i, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j), y(i, j)
                        endif
                        if (ijct(i,j)==4) Then
                            write (27, "(2F12.3, '  -1')") x(i, j), y(i, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j), y(i+1, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j+1), y(i+1, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j), y(i, j)
                        endif
                        if (ijct(i,j)==3) Then
                            write (27, "(2F12.3, '  -1')") x(i+1, j), y(i+1, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j+1), y(i+1, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j+1), y(i, j+1)
                            write (27, "(2F12.3, '   1')") x(i+1, j), y(i+1, j)
                        endif
                        if (ijct(i,j)==2) Then
                            write (27, "(2F12.3, '  -1')") x(i, j), y(i, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j+1), y(i+1, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j+1), y(i, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j), y(i, j)
                        endif
                        if (ijct(i,j)==5) Then
                            write (27, "(2F12.3, '  -1')") x(i, j), y(i, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j), y(i+1, j)
                            write (27, "(2F12.3, '   1')") x(i+1, j+1), y(i+1, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j+1), y(i, j+1)
                            write (27, "(2F12.3, '   1')") x(i, j), y(i, j)
                        endif
                        if (ijct(i,j)==1) Then
                            write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==4) Then
                            write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==3) Then
                            write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (25, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==2) Then
                            write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==5) Then
                            write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (25, "('  0', /, 'SEQEND')")
                        endif
                        icompt = i
                        jcompt = j
                        rmin = 0.
                        write (18, "(4I10, 5X, F12.4)") i, j, icompt, jcompt, rmin
                    endif
                enddo
            enddo
        endif
        write (18, "(4I10, 5X, F12.4)") nwcells
        write (15, "(1X, 'NWCELLS=', I10)") nwcells
        write (6, "(1X, 'NWCELLS=', I10)") nwcells
        write (66, "(1X, 'NWCELLS=', I10)") nwcells
        close (18)
        close (30)
        open (39, file='gridext.out', status='unknown')
        do j = jmino, jmaxo
            do i = imino, imaxo
                if (x(i,j)/=0.0 .And. y(i,j)/=0.0) Then
                    write (39, "(2I5, 2X, F10.6, 2X, F10.6)") i, j, x(i, j), y(i, j)
                endif
            enddo
        enddo
        close (39)
        !********** (ntype>=8) else **********
    else
        if (ntype==1 .Or. ntype==3) Then
            read (4, *) dummy
            read (4, *) dummy
            read (4, *) ilt, jlt, x(ilt, jlt), y(ilt, jlt)
            xl = x(ilt, jlt)
            yl = y(ilt, jlt)
            read (4, *) dummy
            read (4, *) dummy
            read (4, *) ift, jft, x(ift, jft), y(ift, jft)
            write (10, "(1X, E12.4, 5X, E12.4)") x(ift, jft), y(ift, jft)
            xf = x(ift, jft)
            yf = y(ift, jft)
            ibp(1) = ift
            jbp(1) = jft
            ksbp(ift, jft) = 1
            ksgi(ift, jft) = 1
            rsgi(ift, jft) = 1.
            dx = xf - xl
            dy = yf - yl
            ang = atan2(dy, dx)
            write (6, "(1X, 'ANG = ', E11.4)") ang
            write (66, "(1X, 'ANG = ', E11.4)") ang
            do n = 2, nbpp - 1
                nn = n + nbpp - 1
                read (4, *) i, j, x(i, j), y(i, j)
                write (10, "(1X, E12.4, 5X, E12.4)") x(i, j), y(i, j)
                xn(i, j) = (x(i,j)-xl)*cos(ang) + (y(i,j)-yl)*sin(ang)
                yn(i, j) = -(x(i,j)-xl)*sin(ang) + (y(i,j)-yl)*cos(ang)
                ibp(n) = i
                jbp(n) = j
                ksbp(i, j) = 1
                ksgi(i, j) = 1
                rsgi(i, j) = 1.
                ii = 2*ift - i
                ibp(nn) = ii
                jbp(nn) = j
                ksbp(ii, j) = 1
                ksgi(ii, j) = 1
                rsgi(ii, j) = 1.
                xn(ii, j) = xn(i, j)
                yn(ii, j) = -yn(i, j)
                x(i, j) = xl + xn(i, j)*cos(ang) - yn(i, j)*sin(ang)
                y(i, j) = yl + xn(i, j)*sin(ang) + yn(i, j)*cos(ang)
                x(ii, j) = xl + xn(ii, j)*cos(ang) - yn(ii, j)*sin(ang)
                y(ii, j) = yl + xn(ii, j)*sin(ang) + yn(ii, j)*cos(ang)
            enddo
            read (4, *) i, j, x(i, j), y(i, j)
            write (10, "(1X, E12.4, 5X, E12.4)") x(i, j), y(i, j)
            write (10, "(1X, E12.4, 5X, E12.4)") xf, yf
            ibp(nbpp) = i
            jbp(nbpp) = j
            ksbp(i, j) = 1
            ksgi(i, j) = 1
            rsgi(i, j) = 1.
        endif
        if (ntype==2 .Or. ntype==4) Then
            read (4, *) dummy
            read (4, *) dummy
            read (4, *) ilt, jlt, x(ilt, jlt), y(ilt, jlt)
            xl = x(ilt, jlt)
            yl = y(ilt, jlt)
            read (4, *) dummy
            read (4, *) dummy
            read (4, *) ift, jft, x(ift, jft), y(ift, jft)
            write (10, "(1X, E12.4, 5X, E12.4)") x(ift, jft), y(ift, jft)
            xf = x(ift, jft)
            yf = y(ift, jft)
            ibp(1) = ift
            jbp(1) = jft
            ksbp(ift, jft) = 1
            ksgi(ift, jft) = 1
            rsgi(ift, jft) = 1.
            dx = xf - xl
            dy = yf - yl
            ang = atan2(dy, dx)
            do n = 2, nbpp - 1
                nn = n + nbpp - 1
                read (4, *) i, j, x(i, j), y(i, j)
                write (10, "(1X, E12.4, 5X, E12.4)") x(i, j), y(i, j)
                xn(i, j) = (x(i,j)-xl)*cos(ang) + (y(i,j)-yl)*sin(ang)
                yn(i, j) = -(x(i,j)-xl)*sin(ang) + (y(i,j)-yl)*cos(ang)
                ibp(n) = i
                jbp(n) = j
                ksbp(i, j) = 1
                ksgi(i, j) = 1
                rsgi(i, j) = 1.
                jj = 2*jft - j
                ibp(nn) = i
                jbp(nn) = jj
                ksbp(i, jj) = 1
                ksgi(i, jj) = 1
                rsgi(i, jj) = 1.
                xn(i, jj) = xn(i, j)
                yn(i, jj) = -yn(i, j)
                x(i, j) = xl + xn(i, j)*cos(ang) - yn(i, j)*sin(ang)
                y(i, j) = yl + xn(i, j)*sin(ang) + yn(i, j)*cos(ang)
                x(i, jj) = xl + xn(i, jj)*cos(ang) - yn(i, jj)*sin(ang)
                y(i, jj) = yl + xn(i, jj)*sin(ang) + yn(i, jj)*cos(ang)
            enddo
            read (4, *) i, j, x(i, j), y(i, j)
            write (10, "(1X, E12.4, 5X, E12.4)") x(i, j), y(i, j)
            write (10, "(1X, E12.4, 5X, E12.4)") xf, yf
            ibp(nbpp) = i
            jbp(nbpp) = j
            ksbp(i, j) = 1
            ksgi(i, j) = 1
            rsgi(i, j) = 1.
        endif
        if (ntype>=5 .And. ntype<=6) Then
            read (4, *) dummy
            read (4, *) dummy
            read (4, *) ilt, jlt, x(ilt, jlt), y(ilt, jlt)
            write (7, *) ilt, jlt, x(ilt, jlt), y(ilt, jlt)
            xl = x(ilt, jlt)
            yl = y(ilt, jlt)
            read (4, *) dummy
            read (4, *) dummy
            read (4, *) ift, jft, x(ift, jft), y(ift, jft)
            write (7, *) ift, jft, x(ift, jft), y(ift, jft)
            write (10, "(1X, E12.4, 5X, E12.4)") x(ift, jft), y(ift, jft)
            xf = x(ift, jft)
            yf = y(ift, jft)
            xn(ift, jft) = 0.
            yn(ift, jft) = 0.
            ibp(1) = ift
            jbp(1) = jft
            ksbp(ift, jft) = 1
            ksgi(ift, jft) = 1
            rsgi(ift, jft) = 1.
            do n = 2, nbpp - 1
                read (4, *) i, j, x(i, j), y(i, j)
                write (7, *) i, j, x(i, j), y(i, j)
                write (10, "(1X, E12.4, 5X, E12.4)") x(i, j), y(i, j)
                ibp(n) = i
                jbp(n) = j
                ksbp(i, j) = 1
                ksgi(i, j) = 1
                rsgi(i, j) = 1.
            enddo
            read (4, *) i, j, x(i, j), y(i, j)
            write (7, *) i, j, x(i, j), y(i, j)
            write (10, "(1X, E12.4, 5X, E12.4)") x(i, j), y(i, j)
            write (10, "(1X, E12.4, 5X, E12.4)") xf, yf
            ibp(nbpp) = i
            jbp(nbpp) = j
            ksbp(i, j) = 1
            ksgi(i, j) = 1
            rsgi(i, j) = 1.
        endif
        do j = jmin, jmax
            do i = imin, imax
                x(i, j) = x(i, j) + xshift
                y(i, j) = y(i, j) + yshift
            enddo
        enddo
        do is = 1, 4
            read (3, "(130X)")
        enddo
        read (3, "(I3, 2X, 125I1)") jctmp
        close (3)
        open (3, file='cell.inp', status='unknown')
        do is = 1, 4
            read (3, "(130X)")
        enddo
        if (jctmp/=jc) Then
            do jt = 1, jc, 125
                jf = jt
                jlast = jt + 124
                if (jlast>jc) jlast = jc
                write (7, "(1X, '  CELL TYPE ARRAY,J=', I5, 2X, 'TO J=', I5, //)") jf, jlast
                do i = 1, ic
                    read (3, "(120I1)")(ijct(i,j), j=jf, jlast)
                    write (7, "(1X, I3, 2X, 125I1)") i, (ijct(i,j), j=jf, jlast)
                enddo
            enddo
        else
            do it = 1, ic, 125
                ifirst = it
                ilast = it + 124
                if (ilast>ic) ilast = ic
                write (7, "(1X, '  CELL TYPE ARRAY,I=', I5, 2X, 'TO I=', I5, //)") ifirst, ilast
                do j = jc, 1, -1
                    read (3, "(I3, 2X, 125I1)") jdumy, (ijct(i,j), i=ifirst, ilast)
                    write (7, "(1X, I3, 2X, 125I1)") jdumy, (ijct(i,j), i=ifirst, ilast)
                enddo
            enddo
        endif


        !****** (ntype/=0) if ******
        if (ntype/=0) Then
            !********** (ntype/=7) if **********
            if (ntype/=7) then
                do j = jmino + 1, jmaxo
                    do i = imino + 1, imaxo
                        if (ksgi(i,j)==0) Then
                            if (ijct(i,j)>=1 .And. ijct(i,j)<=5) Then
                                ksgi(i, j) = 1
                                rsgi(i, j) = 1.
                            endif
                            if (ijct(i-1,j)>=1 .And. ijct(i-1,j)<=5) Then
                                ksgi(i, j) = 1
                                rsgi(i, j) = 1.
                            endif
                            if (ijct(i,j-1)>=1 .And. ijct(i,j-1)<=5) Then
                                ksgi(i, j) = 1
                                rsgi(i, j) = 1.
                            endif
                        endif
                    enddo
                enddo
                if (ntype==1) Then
                    do i = ift - 1, imin, -1
                        ir = 2*ift - i
                        do j = jmin, jmax
                            ksgi(ir, j) = ksgi(i, j)
                            rsgi(ir, j) = rsgi(i, j)
                        enddo
                    enddo
                endif
                if (ntype==3) Then
                    do i = ift + 1, imax
                        ir = 2*ift - i
                        do j = jmin, jmax
                            ksgi(ir, j) = ksgi(i, j)
                            rsgi(ir, j) = rsgi(i, j)
                        enddo
                    enddo
                endif
                if (ntype==2) Then
                    do i = imin, imax
                        do j = jft - 1, jmin, -1
                            jr = 2*jft - j
                            ksgi(i, jr) = ksgi(i, j)
                            rsgi(i, jr) = rsgi(i, j)
                        enddo
                    enddo
                endif
                if (ntype==4) Then
                    do i = imin, imax
                        do j = jft + 1, jmax
                            jr = 2*jft - j
                            ksgi(i, jr) = ksgi(i, j)
                            rsgi(i, jr) = rsgi(i, j)
                        enddo
                    enddo
                endif
                write (7, "(1X, 'KSGI ARRAY', //)")
                do jt = jmin, jmax, 120
                    jf = jt
                    jlast = jt + 119
                    if (jlast>jmax) jlast = jmax
                    write (7, "(1X, '  CELL TYPE ARRAY,J=', I5, 2X, 'TO J=', I5, //)") jf, jlast
                    do i = imin, imax
                        write (7, "(1X, I3, 2X, 120I1)") i, (ksgi(i,j), j=jf, jlast)
                    enddo
                enddo
                nred = 0
                nblk = 0
                do i = imin, imax
                    do j = jmin, jmax
                        if (ksbp(i,j)==0 .And. ksgi(i,j)==1) Then
                            ipj = i + j
                            ir = mod(ipj, 2)
                            if (ir==0) Then
                                nred = nred + 1
                                ired(nred) = i
                                jred(nred) = j
                            else
                                nblk = nblk + 1
                                iblk(nblk) = i
                                jblk(nblk) = j
                            endif
                        endif
                    enddo
                enddo

                !��X��Y��ʼ��Ϊͨ����ϵ��������˹�����ɳڻ�õ��ڲ�ֵ��RKII��RKJJ�ѱ���ʼ��Ϊ1.0��
                do j = jmin, jmax
                    do i = imin, imax
                        if(i==imin)then
                            if(j==jmin)then
                                cij(i, j) = -(rkii(i,j)+rkjj(i,j))
                            else
                                cij(i, j) = -(rkii(i,j)+rkjj(i,j)+rkjj(i,j-1))
                            endif
                        elseif(j==jmin)then
                            cij(i, j) = -(rkii(i,j)+rkii(i-1,j)+rkjj(i,j))
                        else
                            cij(i,j) = -(rkii(i,j)+rkii(i-1,j)+rkjj(i,j)+rkjj(i,j-1))
                        endif
                    enddo
                enddo

                itr = 0
                do while (.true.)
                    itr = itr + 1
                    rsqx = 0.
                    rsqy = 0.
                    do n = 1, nred
                        i = ired(n)
                        j = jred(n)
                        rsd = rkii(i, j)*x(i+1, j) + rkii(i-1, j)*x(i-1, j) + rkjj(i, j)*x(i, j+1) + rkjj(i, j-1)*x(i, j-1) + cij(i, j)*x(i, j)
                        x(i, j) = x(i, j) - rpx*rsd/cij(i, j)
                        rsqx = rsqx + rsd*rsd
                        rsd = rkii(i, j)*y(i+1, j) + rkii(i-1, j)*y(i-1, j) + rkjj(i, j)*y(i, j+1) + rkjj(i, j-1)*y(i, j-1) + cij(i, j)*y(i, j)
                        y(i, j) = y(i, j) - rpx*rsd/cij(i, j)
                        rsqy = rsqy + rsd*rsd
                    enddo
                    do n = 1, nblk
                        i = iblk(n)
                        j = jblk(n)
                        rsd = rkii(i, j)*x(i+1, j) + rkii(i-1, j)*x(i-1, j) + rkjj(i, j)*x(i, j+1) + rkjj(i, j-1)*x(i, j-1) + cij(i, j)*x(i, j)
                        x(i, j) = x(i, j) - rpx*rsd/cij(i, j)
                        rsqx = rsqx + rsd*rsd
                        rsd = rkii(i, j)*y(i+1, j) + rkii(i-1, j)*y(i-1, j) + rkjj(i, j)*y(i, j+1) + rkjj(i, j-1)*y(i, j-1) + cij(i, j)*y(i, j)
                        y(i, j) = y(i, j) - rpx*rsd/cij(i, j)
                        rsqy = rsqy + rsd*rsd
                    enddo
                    if(itr>=itrxm)exit
                    if (rsqx<rsqxm .and. rsqy<rsqxm)exit
                enddo

                write (6, "(1X, 'DifF INITIAL X&Y, ITER = ', I3, ' RSX,RSY =', 2(2X,E11.4))") itr, rsqx, rsqy
                write (66, "(1X, 'DifF INITIAL X&Y, ITER = ', I3, ' RSX,RSY =', 2(2X,E11.4))") itr, rsqx, rsqy
                if (ntype==1) Then
                    do j = jmin, jmax
                        write (9, "(1X, 'J=', I5)") j
                        write (13, "(1X, 'J=', I5)") j
                        do i = imin, imaxo
                            if (ksgi(i,j)==1) Then
                                write (9, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (13, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                    do i = imin, imaxo
                        write (9, "(1X, 'I=', I5)") i
                        write (13, "(1X, 'I=', I5)") i
                        do j = jmin, jmax
                            if (ksgi(i,j)==1) Then
                                write (8, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (13, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                endif
                if (ntype==2) Then
                    do j = jmin, jmaxo
                        write (9, "(1X, 'J=', I5)") j
                        write (13, "(1X, 'J=', I5)") j
                        do i = imin, imax
                            if (ksgi(i,j)==1) Then
                                write (9, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (13, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                    do i = imin, imax
                        write (9, "(1X, 'I=', I5)") i
                        write (13, "(1X, 'I=', I5)") i
                        do j = jmin, jmaxo
                            if (ksgi(i,j)==1) Then
                                write (8, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (13, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                endif
                if (ntype==3) Then
                    do j = jmin, jmax
                        write (9, "(1X, 'J=', I5)") j
                        write (13, "(1X, 'J=', I5)") j
                        do i = imino, imax
                            if (ksgi(i,j)==1) Then
                                write (9, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (13, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                    do i = imino, imax
                        write (9, "(1X, 'I=', I5)") i
                        write (13, "(1X, 'I=', I5)") i
                        do j = jmin, jmax
                            if (ksgi(i,j)==1) Then
                                write (8, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (13, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                endif
                if (ntype==4) Then
                    do j = jmino, jmax
                        write (9, "(1X, 'J=', I5)") j
                        write (13, "(1X, 'J=', I5)") j
                        do i = imin, imax
                            if (ksgi(i,j)==1) Then
                                write (8, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (13, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                    do i = imin, imax
                        write (9, "(1X, 'I=', I5)") i
                        write (13, "(1X, 'I=', I5)") i
                        do j = jmino, jmax
                            if (ksgi(i,j)==1) Then
                                write (8, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (13, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                endif
                if (ntype>=5) Then
                    do j = jmin, jmax
                        write (9, "(1X, 'J=', I5)") j
                        write (13, "(1X, 'J=', I5)") j
                        do i = imin, imax
                            if (ksgi(i,j)==1) Then
                                write (9, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (13, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                    do i = imin, imax
                        write (9, "(1X, 'I=', I5)") i
                        write (13, "(1X, 'I=', I5)") i
                        do j = jmin, jmax
                            if (ksgi(i,j)==1) Then
                                write (8, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (13, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                endif
                write (7, "(1X, 'INITIAL X AND Y FIELDS', /)")
                write (7, "(5X, 'I', 10X, 'J', 11X, 'X', 12X, 'Y', /)")
                do j = jmin, jmax
                    do i = imin, imax
                        if (ksgi(i,j)==1) Then
                            write (7, "(1X, I10, 5X, I10, 5X, F12.4, 5X, F12.4)") i, j, x(i, j), y(i, j)
                        endif
                    enddo
                enddo
                write (7, "(//////)")
                do j = jmin, jmax
                    do i = imin, imax
                        if (ijct(i,j)==1) Then
                            write (26, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==4) Then
                            write (26, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==3) Then
                            write (26, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (26, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==2) Then
                            write (26, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==5) Then
                            write (26, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'SEQEND')")
                        endif
                    enddo
                enddo
                
                itrg = 0
                !********** (ntype/=6) if **********
                if (ntype/=6) then
                    do while(.true.)
                        itrg = itrg + 1
                        do n = 1, nbp
                            i = ibp(n)
                            j = jbp(n)
                            if (ksgi(i-1,j)==0 .And. ksgi(i+1,j)==1) Then
                                if (ksgi(i+2,j)==0) Then
                                    dxdi = x(i+1, j) - x(i, j)
                                    dydi = y(i+1, j) - y(i, j)
                                else
                                    dxdi = 0.5*(4.*x(i+1,j)-3.*x(i,j)-x(i+2,j))
                                    dydi = 0.5*(4.*y(i+1,j)-3.*y(i,j)-y(i+2,j))
                                endif
                            endif
                            if (ksgi(i+1,j)==0 .And. ksgi(i-1,j)==1) Then
                                if (ksgi(i-2,j)==0) Then
                                    dxdi = x(i, j) - x(i-1, j)
                                    dydi = y(i, j) - y(i-1, j)
                                else
                                    dxdi = 0.5*(3.*x(i,j)-4.*x(i-1,j)+x(i-2,j))
                                    dydi = 0.5*(3.*y(i,j)-4.*y(i-1,j)+y(i-2,j))
                                endif
                            endif
                            if (ksgi(i+1,j)==1 .And. ksgi(i-1,j)==1) Then
                                dxdi = 0.5*(x(i+1,j)-x(i-1,j))
                                dydi = 0.5*(y(i+1,j)-y(i-1,j))
                            endif
                            if (ksgi(i,j-1)==0 .And. ksgi(i,j+1)==1) Then
                                if (ksgi(i,j+2)==0) Then
                                    dxdj = x(i, j+1) - x(i, j)
                                    dydj = y(i, j+1) - y(i, j)
                                else
                                    dxdj = 0.5*(4.*x(i,j+1)-3.*x(i,j)-x(i,j+2))
                                    dydj = 0.5*(4.*y(i,j+1)-3.*y(i,j)-y(i,j+2))
                                endif
                            endif
                            if (ksgi(i,j+1)==0 .And. ksgi(i,j-1)==1) Then
                                if (ksgi(i,j-2)==0) Then
                                    dxdj = x(i, j) - x(i, j-1)
                                    dydj = y(i, j) - y(i, j-1)
                                else
                                    dxdj = 0.5*(3.*x(i,j)-4.*x(i,j-1)+x(i,j-2))
                                    dydj = 0.5*(3.*y(i,j)-4.*y(i,j-1)+y(i,j-2))
                                endif
                            endif
                            if (ksgi(i,j+1)==1 .And. ksgi(i,j-1)==1) Then
                                dxdj = 0.5*(x(i,j+1)-x(i,j-1))
                                dydj = 0.5*(y(i,j+1)-y(i,j-1))
                            endif
                            hi(i, j) = sqrt(dxdi*dxdi+dydi*dydi)
                            hj(i, j) = sqrt(dxdj*dxdj+dydj*dydj)
                            rki(i, j) = hj(i, j)/hi(i, j)
                        enddo
                        if (isirki==1) Then
                            if (jsirki==1) Then
                                do j = jmin, jmax
                                    do i = imin, imax
                                        rkii(i, j) = 1.
                                        rkjj(i, j) = rkjdki
                                    enddo
                                enddo
                            endif
                            itr = 0
                            do j = jmin, jmax
                                do i = imin, imax
                                    if(i==imin)then
                                        if(j==jmin)then
                                            cij(i, j) = -(rkii(i,j)+rkjj(i,j))
                                        else
                                            cij(i, j) = -(rkii(i,j)+rkjj(i,j)+rkjj(i,j-1))
                                        endif
                                    elseif(j==jmin)then
                                        cij(i, j) = -(rkii(i,j)+rkii(i-1,j)+rkjj(i,j))
                                    else
                                        cij(i,j) = -(rkii(i,j)+rkii(i-1,j)+rkjj(i,j)+rkjj(i,j-1))
                                    endif
                                enddo
                            enddo

                            do while (.true.)
                                itr = itr + 1
                                rsq = 0.
                                do n = 1, nred
                                    i = ired(n)
                                    j = jred(n)
                                    rsd = rkii(i, j)*rki(i+1, j) + rkii(i-1, j)*rki(i-1, j) + rkjj(i, j)*rki(i, j+1) + rkjj(i, j-1)*rki(i, j-1) + cij(i, j)*rki(i, j)
                                    rki(i, j) = rki(i, j) - rpk*rsd/cij(i, j)
                                    rsq = rsq + rsd*rsd
                                enddo
                                do n = 1, nblk
                                    i = iblk(n)
                                    j = jblk(n)
                                    rsd = rkii(i, j)*rki(i+1, j) + rkii(i-1, j)*rki(i-1, j) + rkjj(i, j)*rki(i, j+1) + rkjj(i, j-1)*rki(i, j-1) + cij(i, j)*rki(i, j)
                                    rki(i, j) = rki(i, j) - rpk*rsd/cij(i, j)
                                    rsq = rsq + rsd*rsd
                                enddo
                                if (itr>=itrkm) exit
                                if (rsq<rsqkm) exit
                            enddo
                            write (6, "(1X, 'DifFuse RKI, ITERATION = ', I3, ' RSK =', 2X, E11.4)") itr, rsq
                            write (66, "(1X, 'DifFuse RKI, ITERATION = ', I3, ' RSK =', 2X, E11.4)") itr, rsq
                            do j = jmin, jmax
                                do i = imin, imax
                                    rkj(i, j) = rkjdki/rki(i, j)
                                enddo
                            enddo
                            do j = jmin, jmax
                                do i = imin, imax - 1
                                    rkii(i, j) = 0.5*(rki(i+1,j)+rki(i,j))
                                enddo
                            enddo
                            do j = jmin, jmax - 1
                                do i = imin, imax
                                    rkjj(i, j) = 0.5*(rkj(i,j+1)+rkj(i,j))
                                enddo
                            enddo
                        endif
                        if (isihihj==1) Then
                            if (jsihihj==1) Then
                                do j = jmin, jmax
                                    do i = imin, imax
                                        rkii(i, j) = 1.
                                        rkjj(i, j) = rkjdki
                                    enddo
                                enddo
                            endif
                            do j = jmin, jmax
                                do i = imin, imax
                                    if(i==imin)then
                                        if(j==jmin)then
                                            cij(i, j) = -(rkii(i,j)+rkjj(i,j))
                                        else
                                            cij(i, j) = -(rkii(i,j)+rkjj(i,j)+rkjj(i,j-1))
                                        endif
                                    elseif(j==jmin)then
                                        cij(i, j) = -(rkii(i,j)+rkii(i-1,j)+rkjj(i,j))
                                    else
                                        cij(i,j) = -(rkii(i,j)+rkii(i-1,j)+rkjj(i,j)+rkjj(i,j-1))
                                    endif
                                enddo
                            enddo

                            itr = 0
                            do while (.true.)
                                itr = itr + 1
                                rsqhi = 0.
                                rsqhj = 0.
                                do n = 1, nred
                                    i = ired(n)
                                    j = jred(n)
                                    rsd = rkii(i, j)*hi(i+1, j) + rkii(i-1, j)*hi(i-1, j) + rkjj(i, j)*hi(i, j+1) + rkjj(i, j-1)*hi(i, j-1) + cij(i, j)*hi(i, j)
                                    hi(i, j) = hi(i, j) - rph*rsd/cij(i, j)
                                    rsqhi = rsqhi + rsd*rsd
                                    rsd = rkii(i, j)*hj(i+1, j) + rkii(i-1, j)*hj(i-1, j) + rkjj(i, j)*hj(i, j+1) + rkjj(i, j-1)*hj(i, j-1) + cij(i, j)*hj(i, j)
                                    hj(i, j) = hj(i, j) - rph*rsd/cij(i, j)
                                    rsqhj = rsqhj + rsd*rsd
                                enddo
                                do n = 1, nblk
                                    i = iblk(n)
                                    j = jblk(n)
                                    rsd = rkii(i, j)*hi(i+1, j) + rkii(i-1, j)*hi(i-1, j) + rkjj(i, j)*hi(i, j+1) + rkjj(i, j-1)*hi(i, j-1) + cij(i, j)*hi(i, j)
                                    hi(i, j) = hi(i, j) - rpx*rsd/cij(i, j)
                                    rsqhi = rsqhi + rsd*rsd
                                    rsd = rkii(i, j)*hj(i+1, j) + rkii(i-1, j)*hj(i-1, j) + rkjj(i, j)*hj(i, j+1) + rkjj(i, j-1)*hj(i, j-1) + cij(i, j)*hj(i, j)
                                    hj(i, j) = hj(i, j) - rpx*rsd/cij(i, j)
                                    rsqhj = rsqhj + rsd*rsd
                                enddo
                                if (itr>=itrhm) exit
                                if (rsqhi<rsqhm .and. rsqhj<rsqhm) exit
                            enddo

                            write (6, "(1X, 'DifF HI & HJ, ITER = ', I3, ' RSI,RSJ =', 2(2X,E11.4))") itr, rsqhi, rsqhj
                            write (66, "(1X, 'DifF HI & HJ, ITER = ', I3, ' RSI,RSJ =', 2(2X,E11.4))") itr, rsqhi, rsqhj
                            do j = jmin, jmax
                                do i = imin, imax
                                    rki(i, j) = hj(i, j)/hi(i, j)
                                    rkj(i, j) = hi(i, j)/hj(i, j)
                                enddo
                            enddo
                            do j = jmin, jmax
                                do i = imin, imax - 1
                                    rkii(i, j) = 0.5*(rki(i+1,j)+rki(i,j))
                                enddo
                            enddo
                            do j = jmin, jmax - 1
                                do i = imin, imax
                                    rkjj(i, j) = 0.5*(rkj(i,j+1)+rkj(i,j))
                                enddo
                            enddo
                        endif
                        do j = jmin, jmax
                            do i = imin, imax
                                if(i==imin)then
                                    if(j==jmin)then
                                        cij(i, j) = -(rkii(i,j)+rkjj(i,j))
                                    else
                                        cij(i, j) = -(rkii(i,j)+rkjj(i,j)+rkjj(i,j-1))
                                    endif
                                elseif(j==jmin)then
                                    cij(i, j) = -(rkii(i,j)+rkii(i-1,j)+rkjj(i,j))
                                else
                                    cij(i,j) = -(rkii(i,j)+rkii(i-1,j)+rkjj(i,j)+rkjj(i,j-1))
                                endif
                            enddo
                        enddo

                        itr = 0
                        do while (.true.)
                            itr = itr + 1
                            rsqx = 0.
                            rsqy = 0.
                            do n = 1, nred
                                i = ired(n)
                                j = jred(n)
                                rsd = rkii(i, j)*x(i+1, j) + rkii(i-1, j)*x(i-1, j) + rkjj(i, j)*x(i, j+1) + rkjj(i, j-1)*x(i, j-1) + cij(i, j)*x(i, j)
                                x(i, j) = x(i, j) - rpx*rsd/cij(i, j)
                                rsqx = rsqx + rsd*rsd
                                rsd = rkii(i, j)*y(i+1, j) + rkii(i-1, j)*y(i-1, j) + rkjj(i, j)*y(i, j+1) + rkjj(i, j-1)*y(i, j-1) + cij(i, j)*y(i, j)
                                y(i, j) = y(i, j) - rpx*rsd/cij(i, j)
                                rsqy = rsqy + rsd*rsd
                            enddo
                            do n = 1, nblk
                                i = iblk(n)
                                j = jblk(n)
                                rsd = rkii(i, j)*x(i+1, j) + rkii(i-1, j)*x(i-1, j) + rkjj(i, j)*x(i, j+1) + rkjj(i, j-1)*x(i, j-1) + cij(i, j)*x(i, j)
                                x(i, j) = x(i, j) - rpx*rsd/cij(i, j)
                                rsqx = rsqx + rsd*rsd
                                rsd = rkii(i, j)*y(i+1, j) + rkii(i-1, j)*y(i-1, j) + rkjj(i, j)*y(i, j+1) + rkjj(i, j-1)*y(i, j-1) + cij(i, j)*y(i, j)
                                y(i, j) = y(i, j) - rpx*rsd/cij(i, j)
                                rsqy = rsqy + rsd*rsd
                            enddo
                            if (itr>=itrxm) exit
                            if (rsqx<rsqxm .and. rsqy<rsqxm) exit
                        enddo

                        write (6, "(1X, 'DifF X & Y, ITER = ', I3, ' RSX,RSY =', 2(2X,E11.4))") itr, rsqx, rsqy
                        write (66, "(1X, 'DifF X & Y, ITER = ', I3, ' RSX,RSY =', 2(2X,E11.4))") itr, rsqx, rsqy

                        if (isirki==1) Then
                            ajbm = 0.
                            agijma = 0.
                            agijmi = 1.E+10
                            angmax = 0.
                            angmin = 1.E+10
                            rsqki = 0.
                            do n = 1, nred
                                i = ired(n)
                                j = jred(n)
                                dxdi = 0.5*(x(i+1,j)-x(i-1,j))
                                dydi = 0.5*(y(i+1,j)-y(i-1,j))
                                dxdj = 0.5*(x(i,j+1)-x(i,j-1))
                                dydj = 0.5*(y(i,j+1)-y(i,j-1))
                                gii = dxdi*dxdi + dydi*dydi
                                gjj = dxdj*dxdj + dydj*dydj
                                gij = dxdi*dxdj + dydi*dydj
                                agij = abs(1.-gij*gij/(gii*gjj))
                                angd = (acos(gij/sqrt(gii*gjj)))*57.29578
                                angerr = abs(90.-angd)
                                angmax = max(angmax, angerr)
                                angmin = min(angmin, angerr)
                                hical = sqrt(gii)
                                hjcal = sqrt(gjj)
                                hihjcal = hical*hjcal
                                rkical = hjcal/hical
                                rsqki = rsqki + (rkical-rki(i,j))**2
                            enddo
                            do n = 1, nblk
                                i = iblk(n)
                                j = jblk(n)
                                dxdi = 0.5*(x(i+1,j)-x(i-1,j))
                                dydi = 0.5*(y(i+1,j)-y(i-1,j))
                                dxdj = 0.5*(x(i,j+1)-x(i,j-1))
                                dydj = 0.5*(y(i,j+1)-y(i,j-1))
                                gii = dxdi*dxdi + dydi*dydi
                                gjj = dxdj*dxdj + dydj*dydj
                                gij = dxdi*dxdj + dydi*dydj
                                agij = abs(1.-gij*gij/(gii*gjj))
                                angd = (acos(gij/sqrt(gii*gjj)))*57.29578
                                angerr = abs(90.-angd)
                                angmax = max(angmax, angerr)
                                angmin = min(angmin, angerr)
                                hical = sqrt(gii)
                                hjcal = sqrt(gjj)
                                hihjcal = hical*hjcal
                                rkical = hjcal/hical
                                rsqki = rsqki + (rkical-rki(i,j))**2
                            enddo
                            write (6, "(1X, 'GRID GENERATION LOOP ITERATION =', I10)") itrg
                            write (6, "(1X, 'GLOBAL RES SQ DifF IN RKI=', E12.4)") rsqki
                            write (6, "(1X, 'MIN AND MAX DEVIATION FROM ORTHO =', 2(2X,E12.4))") angmin, angmax
                            write (66, "(1X, 'GRID GENERATION LOOP ITERATION =', I10)") itrg
                            write (66, "(1X, 'GLOBAL RES SQ DifF IN RKI=', E12.4)") rsqki
                            write (66, "(1X, 'MIN AND MAX DEVIATION FROM ORTHO =', 2(2X,E12.4))") angmin, angmax
                            if (itrg>=itrgm) exit
                            if (angmax>=angoro) cycle
                        endif

                        if (isihihj==1) Then
                            rsqhi = 0.
                            rsqhj = 0.
                            ajbm = 0.
                            agijma = 0.
                            agijmi = 1.E+10
                            angmax = 0.
                            angmin = 1.E+10
                            do n = 1, nred
                                i = ired(n)
                                j = jred(n)
                                dxdi = 0.5*(x(i+1,j)-x(i-1,j))
                                dydi = 0.5*(y(i+1,j)-y(i-1,j))
                                dxdj = 0.5*(x(i,j+1)-x(i,j-1))
                                dydj = 0.5*(y(i,j+1)-y(i,j-1))
                                gii = dxdi*dxdi + dydi*dydi
                                gjj = dxdj*dxdj + dydj*dydj
                                gij = dxdi*dxdj + dydi*dydj
                                agij = abs(1.-gij*gij/(gii*gjj))
                                angd = (acos(gij/sqrt(gii*gjj)))*57.29578
                                angerr = abs(90.-angd)
                                angmax = max(angmax, angerr)
                                angmin = min(angmin, angerr)
                                hical = sqrt(gii)
                                hjcal = sqrt(gjj)
                                rsqhi = rsqhi + (hical-hi(i,j))**2
                                rsqhj = rsqhj + (hjcal-hj(i,j))**2
                            enddo
                            do n = 1, nblk
                                i = iblk(n)
                                j = jblk(n)
                                dxdi = 0.5*(x(i+1,j)-x(i-1,j))
                                dydi = 0.5*(y(i+1,j)-y(i-1,j))
                                dxdj = 0.5*(x(i,j+1)-x(i,j-1))
                                dydj = 0.5*(y(i,j+1)-y(i,j-1))
                                gii = dxdi*dxdi + dydi*dydi
                                gjj = dxdj*dxdj + dydj*dydj
                                gij = dxdi*dxdj + dydi*dydj
                                agij = abs(1.-gij*gij/(gii*gjj))
                                angd = (acos(gij/sqrt(gii*gjj)))*57.29578
                                angerr = abs(90.-angd)
                                angmax = max(angmax, angerr)
                                angmin = min(angmin, angerr)
                                hical = sqrt(gii)
                                hjcal = sqrt(gjj)
                                rsqhi = rsqhi + (hical-hi(i,j))**2
                                rsqhj = rsqhj + (hjcal-hj(i,j))**2
                            enddo
                            write (6, "(1X, 'GRID GENERATION LOOP ITERATION =', I10)") itrg
                            write (6, "(1X, 'GLOBAL RES SQ DifF IN HI,HJ=', 2(2X,E12.4))") rsqhi, rsqhj
                            write (6, "(1X, 'MIN AND MAX DEVIATION FROM ORTHO =', 2(2X,E12.4))") angmin, angmax
                            write (66, "(1X, 'GRID GENERATION LOOP ITERATION =', I10)") itrg
                            write (66, "(1X, 'GLOBAL RES SQ DifF IN HI,HJ=', 2(2X,E12.4))") rsqhi, rsqhj
                            write (66, "(1X, 'MIN AND MAX DEVIATION FROM ORTHO =', 2(2X,E12.4))") angmin, angmax
                            if (itrg>=itrgm) exit
                            if (angmax<angoro) exit
                        endif
                    enddo
                else
                !********** (ntype/=6) else **********
                    do while (.true.)
                        itrg = itrg + 1
                        do j = jmin, jmax
                            do i = imin, imax
                                if(i==imin)then
                                    if(j==jmin)then
                                        cij(i, j) = -(rkii(i,j)+rkjj(i,j))
                                    else
                                        cij(i, j) = -(rkii(i,j)+rkjj(i,j)+rkjj(i,j-1))
                                    endif
                                elseif(j==jmin)then
                                    cij(i, j) = -(rkii(i,j)+rkii(i-1,j)+rkjj(i,j))
                                else
                                    cij(i,j) = -(rkii(i,j)+rkii(i-1,j)+rkjj(i,j)+rkjj(i,j-1))
                                endif
                            enddo
                        enddo

                        itr = 0
                        do while (.true.)
                            itr = itr + 1
                            rsqx = 0.
                            rsqy = 0.
                            do n = 1, nred
                                i = ired(n)
                                j = jred(n)
                                rsd = rkii(i, j)*x(i+1, j) + rkii(i-1, j)*x(i-1, j) + rkjj(i, j)*x(i, j+1) + rkjj(i, j-1)*x(i, j-1) + cij(i, j)*x(i, j)
                                x(i, j) = x(i, j) - rpx*rsd/cij(i, j)
                                rsqx = rsqx + rsd*rsd
                                rsd = rkii(i, j)*y(i+1, j) + rkii(i-1, j)*y(i-1, j) + rkjj(i, j)*y(i, j+1) + rkjj(i, j-1)*y(i, j-1) + cij(i, j)*y(i, j)
                                y(i, j) = y(i, j) - rpx*rsd/cij(i, j)
                                rsqy = rsqy + rsd*rsd
                            enddo
                            do n = 1, nblk
                                i = iblk(n)
                                j = jblk(n)
                                rsd = rkii(i, j)*x(i+1, j) + rkii(i-1, j)*x(i-1, j) + rkjj(i, j)*x(i, j+1) + rkjj(i, j-1)*x(i, j-1) + cij(i, j)*x(i, j)
                                x(i, j) = x(i, j) - rpx*rsd/cij(i, j)
                                rsqx = rsqx + rsd*rsd
                                rsd = rkii(i, j)*y(i+1, j) + rkii(i-1, j)*y(i-1, j) + rkjj(i, j)*y(i, j+1) + rkjj(i, j-1)*y(i, j-1) + cij(i, j)*y(i, j)
                                y(i, j) = y(i, j) - rpx*rsd/cij(i, j)
                                rsqy = rsqy + rsd*rsd
                            enddo
                            if (itr>=itrxm) exit
                            if (rsqx<rsqxm .and. rsqy<rsqxm) exit
                        enddo

                        write (6, "(1X, 'DifF X & Y, ITER = ', I3, ' RSX,RSY =', 2(2X,E11.4))") itr, rsqx, rsqy
                        write (66, "(1X, 'DifF X & Y, ITER = ', I3, ' RSX,RSY =', 2(2X,E11.4))") itr, rsqx, rsqy
                        agijma = 0.
                        agijmi = 1.E+10
                        angmax = 0.
                        angmin = 1.E+10
                        do j = jmino, jmaxo - 1
                            do i = imino, imaxo - 1
                                if (ijct(i,j)>=1 .And. ijct(i,j)<=5) Then
                                    dxdi = 0.5*(x(i+1,j+1)-x(i,j+1)+x(i+1,j)-x(i,j))
                                    dydi = 0.5*(y(i+1,j+1)-y(i,j+1)+y(i+1,j)-y(i,j))
                                    dxdj = 0.5*(x(i+1,j+1)-x(i+1,j)+x(i,j+1)-x(i,j))
                                    dydj = 0.5*(y(i+1,j+1)-y(i+1,j)+y(i,j+1)-y(i,j))
                                    gii = dxdi*dxdi + dydi*dydi
                                    gjj = dxdj*dxdj + dydj*dydj
                                    gij = dxdi*dxdj + dydi*dydj
                                    rki(i, j) = gjj
                                    rkj(i, j) = gii
                                    angd = (acos(gij/sqrt(gii*gjj)))*57.29578
                                    angerr = abs(90.-angd)
                                    angmax = max(angmax, angerr)
                                    angmin = min(angmin, angerr)
                                endif
                            enddo
                        enddo
                        do j = jmin + 1, jmax - 1
                            do i = imin, imax - 1
                                rkii(i, j) = 0.5*(rki(i,j)+rki(i,j-1))
                            enddo
                        enddo
                        do j = jmin, jmax - 1
                            do i = imin + 1, imax
                                rkjj(i, j) = 0.5*(rkj(i,j)+rkj(i-1,j))
                            enddo
                        enddo
                        write (6, "(1X, 'GRID GENERATION LOOP ITERATION =', I10)") itrg
                        write (6, "(1X, 'MIN AND MAX DEVIATION FROM ORTHO =', 2(2X,E12.4))") angmin, angmax
                        write (66, "(1X, 'GRID GENERATION LOOP ITERATION =', I10)") itrg
                        write (66, "(1X, 'MIN AND MAX DEVIATION FROM ORTHO =', 2(2X,E12.4))") angmin, angmax
                        if (itrg>=itrgm) exit
                        if (angmax<angoro) exit
                    enddo
                endif
                !********** (ntype/=6) endif **********
            else
            !********** (ntype/=7) elseif **********
                itn7 = 0
                x(ib, jb) = xibjb
                y(ib, jb) = yibjb
                x(ie, jb) = xiejb
                y(ie, jb) = yiejb
                x(ib, je) = xibje
                y(ib, je) = yibje
                x(ie, je) = xieje
                y(ie, je) = yieje
                xn(ib, jb) = xibjb
                yn(ib, jb) = yibjb
                xn(ie, jb) = xiejb
                yn(ie, jb) = yiejb
                xn(ib, je) = xibje
                yn(ib, je) = yibje
                xn(ie, je) = xieje
                yn(ie, je) = yieje
                dytmp = (yibje-yibjb)/float(je-jb)
                ytmp = yibjb
                do j = jb + 1, je - 1
                    ytmp = ytmp + dytmp
                    y(ib, j) = ytmp
                    x(ib, j) = fib(ytmp, j)
                enddo
                dytmp = (yieje-yiejb)/float(je-jb)
                ytmp = yiejb
                do j = jb + 1, je - 1
                    ytmp = ytmp + dytmp
                    y(ie, j) = ytmp
                    x(ie, j) = fie(ytmp, j)
                enddo
                dxtmp = (xiejb-xibjb)/float(ie-ib)
                xtmp = xibjb
                do i = ib + 1, ie - 1
                    xtmp = xtmp + dxtmp
                    x(i, jb) = xtmp
                    y(i, jb) = gjb(xtmp, i)
                enddo
                dxtmp = (xieje-xibje)/float(ie-ib)
                xtmp = xibje
                do i = ib + 1, ie - 1
                    xtmp = xtmp + dxtmp
                    x(i, je) = xtmp
                    y(i, je) = gje(xtmp, i)
                enddo
                open (90, file='ibndry.out', status='unknown')
                write (90, "(' ALONG J=JB ', /)")
                j = jb
                do i = ib, ie
                    write (90, "(2I5, 2F12.4)") i, j, x(i, j), y(i, j)
                enddo
                write (90, "(/)")
                write (90, "(' ALONG I=IE ', /)")
                i = ie
                do j = jb, je
                    write (90, "(2I5, 2F12.4)") i, j, x(i, j), y(i, j)
                enddo
                write (90, "(/)")
                write (90, "(' ALONG J=JE ', /)")
                j = je
                do i = ie, ib, -1
                    write (90, "(2I5, 2F12.4)") i, j, x(i, j), y(i, j)
                enddo
                write (90, "(/)")
                write (90, "(' ALONG I=IB ', /)")
                i = ib
                do j = je, jb, -1
                    write (90, "(2I5, 2F12.4)") i, j, x(i, j), y(i, j)
                enddo
                close (90)
                smd = 1.
                ssq = smd*smd
                cr7 = .5/(1.+ssq)
                nr = 0
                do while (.true.)
                    nr = nr + 1
                    do j = jb + 1, je - 1
                        do i = ib + 1, ie - 1
                            xn(i, j) = cr7*(x(i+1,j)+x(i-1,j)+ssq*(x(i,j+1)+x(i,j-1)))
                            yn(i, j) = cr7*(y(i+1,j)+y(i-1,j)+ssq*(y(i,j+1)+y(i,j-1)))
                        enddo
                    enddo
                    do j = jb + 1, je - 1
                        do i = ib + 1, ie - 1
                            x(i, j) = xn(i, j)
                            y(i, j) = yn(i, j)
                        enddo
                    enddo
                    if (nr>=n7relax) exit
                enddo
                write (6, "(' COMPLETED INIT, NR  = ', I10/)") nr
                Call calsmd
                smdold = smd
                ssq = smd*smd
                smdi = 1./smd
                cr7 = .5/(1.+ssq)
                write (6, "(' SMdoLD,SMD AFTER INIT     = ', 2E14.5/)") smdold, smd
                do j = 1, je
                    do i = 1, ie
                        if (ijct(i,j)==1) Then
                            write (26, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==4) Then
                            write (26, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==3) Then
                            write (26, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (26, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==2) Then
                            write (26, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'SEQEND')")
                        endif
                        if (ijct(i,j)==5) Then
                            write (26, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                            write (26, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                            write (26, "('  0', /, 'SEQEND')")
                        endif
                    enddo
                enddo
                do while (.true.)
                    itn7 = itn7 + 1
                    do j = jb + 1, je - 1
                        x(ib, j) = fib(y(ib,j), j)
                        x(ie, j) = fie(y(ie,j), j)
                    enddo
                    do nxy = 1, nxyit
                        do i = ib + 1, ie - 1
                            if (mod(i+jb,2)==0) Then
                                resd7 = x(i, jb) - x(i, jb+1) - 0.25*smdi*(y(i+1,jb)-y(i-1,jb)+y(i+1,jb+1)-y(i-1,jb+1))
                                x(i, jb) = x(i, jb) - rp7*resd7
                            endif
                            if (mod(i+je,2)==0) Then
                                resd7 = x(i, je) - x(i, je-1) + 0.25*smdi*(y(i+1,je)-y(i-1,je)+y(i+1,je-1)-y(i-1,je-1))
                                x(i, je) = x(i, je) - rp7*resd7
                            endif
                        enddo
                        do j = jb + 1, je - 1
                            do i = ib + 1, ie - 1
                                if (mod(i+j,2)==0) Then
                                    resd7 = x(i, j) - cr7*(x(i+1,j)+x(i-1,j)+ssq*(x(i,j+1)+x(i,j-1)))
                                    x(i, j) = x(i, j) - rp7*resd7
                                endif
                            enddo
                        enddo
                        do i = ib + 1, ie - 1
                            if (mod(i+jb,2)/=0) Then
                                resd7 = x(i, jb) - x(i, jb+1) - 0.25*smdi*(y(i+1,jb)-y(i-1,jb)+y(i+1,jb+1)-y(i-1,jb+1))
                                x(i, jb) = x(i, jb) - rp7*resd7
                            endif
                            if (mod(i+je,2)/=0) Then
                                resd7 = x(i, je) - x(i, je-1) + 0.25*smdi*(y(i+1,je)-y(i-1,je)+y(i+1,je-1)-y(i-1,je-1))
                                x(i, je) = x(i, je) - rp7*resd7
                            endif
                        enddo
                        do j = jb + 1, je - 1
                            do i = ib + 1, ie - 1
                                if (mod(i+j,2)/=0) Then
                                    resd7 = x(i, j) - cr7*(x(i+1,j)+x(i-1,j)+ssq*(x(i,j+1)+x(i,j-1)))
                                    x(i, j) = x(i, j) - rp7*resd7
                                endif
                            enddo
                        enddo
                    enddo
                    smdold = smd
                    Call calsmd
                    ssq = smd*smd
                    smdi = 1./smd
                    cr7 = .5/(1.+ssq)
                    write (6, "(' SMdoLD,SMD AFTER X SWEEP  = ', 2E14.5/)") smdold, smd
                    do i = ib + 1, ie - 1
                        y(i, jb) = gjb(x(i,jb), i)
                        y(i, je) = gje(x(i,je), i)
                    enddo
                    do nxy = 1, nxyit
                        do j = jb + 1, je - 1
                            if (mod(ib+j,2)==0) Then
                                resd7 = y(ib, j) - y(ib+1, j) - 0.25*smd*(x(ib,j+1)-x(ib,j-1)+x(ib+1,j+1)-x(ib+1,j-1))
                                y(ib, j) = y(ib, j) - rp7*resd7
                            endif
                            if (mod(ie+j,2)==0) Then
                                resd7 = y(ie, j) - y(ie-1, j) + 0.25*smd*(x(ie,j+1)-x(ie,j-1)+x(ie-1,j+1)-x(ie-1,j-1))
                                y(ie, j) = y(ie, j) - rp7*resd7
                            endif
                        enddo
                        do j = jb + 1, je - 1
                            do i = ib + 1, ie - 1
                                if (mod(i+j,2)==0) Then
                                    resd7 = y(i, j) - cr7*(y(i+1,j)+y(i-1,j)+ssq*(y(i,j+1)+y(i,j-1)))
                                    y(i, j) = y(i, j) - rp7*resd7
                                endif
                            enddo
                        enddo
                        do j = jb + 1, je - 1
                            if (mod(ib+j,2)/=0) Then
                                resd7 = y(ib, j) - y(ib+1, j) - 0.25*smd*(x(ib,j+1)-x(ib,j-1)+x(ib+1,j+1)-x(ib+1,j-1))
                                y(ib, j) = y(ib, j) - rp7*resd7
                            endif
                            if (mod(ie+j,2)/=0) Then
                                resd7 = y(ie, j) - y(ie-1, j) + 0.25*smd*(x(ie,j+1)-x(ie,j-1)+x(ie-1,j+1)-x(ie-1,j-1))
                                y(ie, j) = y(ie, j) - rp7*resd7
                            endif
                        enddo
                        do j = jb + 1, je - 1
                            do i = ib + 1, ie - 1
                                if (mod(i+j,2)/=0) Then
                                    resd7 = y(i, j) - cr7*(y(i+1,j)+y(i-1,j)+ssq*(y(i,j+1)+y(i,j-1)))
                                    y(i, j) = y(i, j) - rp7*resd7
                                endif
                            enddo
                        enddo
                    enddo
                    smdold = smd
                    Call calsmd
                    ssq = smd*smd
                    smdi = 1./smd
                    cr7 = .5/(1.+ssq)
                    write (6, "(' SMdoLD,SMD AFTER Y SWEEP  = ', 2E14.5/)") smdold, smd
                    serr = abs((smdold-smd)/smd)
                    agijma = 0.
                    agijmi = 1.E+10
                    angmax = 0.
                    angmin = 1.E+10
                    do j = jb, je - 1
                        do i = ib, ie - 1
                            dxdi = 0.5*(x(i+1,j+1)-x(i,j+1)+x(i+1,j)-x(i,j))
                            dydi = 0.5*(y(i+1,j+1)-y(i,j+1)+y(i+1,j)-y(i,j))
                            dxdj = 0.5*(x(i+1,j+1)-x(i+1,j)+x(i,j+1)-x(i,j))
                            dydj = 0.5*(y(i+1,j+1)-y(i+1,j)+y(i,j+1)-y(i,j))
                            gii = dxdi*dxdi + dydi*dydi
                            gjj = dxdj*dxdj + dydj*dydj
                            gij = dxdi*dxdj + dydi*dydj
                            angd = (acos(gij/sqrt(gii*gjj)))*57.29578
                            angerr = abs(90.-angd)
                            angmax = max(angmax, angerr)
                            angmin = min(angmin, angerr)
                        enddo
                    enddo
                    write (6, "(' NTYPE 7, ITER,ITERMAX = ', 2I10/)") itn7, itn7max
                    write (6, "(' SMdoLD,SMD,SERR,SERRMAX = ', 4E14.5/)") smdold, smd, serr, serrmax
                    write (6, "(' ANGMIN,ANGMAX = ', 2E14.4/)") angmin, angmax
                    if (itn7>itn7max) exit
                    if (angmax<=angoro) exit
                    if (serr<=serrmax) exit
                enddo

                x(21, 3) = 13.35
                y(21, 3) = 8.68
                x(22, 3) = 15.14
                y(22, 3) = 6.7
                x(23, 3) = 17.9
                y(23, 3) = 5.8
                x(24, 3) = 20.08
                y(24, 3) = 5.97
                x(22, 2) = 15.1
                y(22, 2) = 3.8
                x(23, 2) = 16.0
                y(23, 2) = 3.5
                open (90, file='shorebndry', status='unknown')
                open (91, file='shoremask', status='unknown')
                i = ib
                do j = jb, je
                    write (90, "(2I5, 2F12.4)") i, j, x(i, j), y(i, j)
                    write (91, "(2F12.3)") x(i, j), y(i, j)
                enddo
                j = je
                do i = ib + 1, ie
                    write (90, "(2I5, 2F12.4)") i, j, x(i, j), y(i, j)
                    write (91, "(2F12.3)") x(i, j), y(i, j)
                enddo
                i = ie
                do j = je - 1, jb, -1
                    write (90, "(2I5, 2F12.4)") i, j, x(i, j), y(i, j)
                    write (91, "(2F12.3)") x(i, j), y(i, j)
                enddo
                j = jb
                do i = ie - 1, 24, -1
                    write (90, "(2I5, 2F12.4)") i, j, x(i, j), y(i, j)
                    write (91, "(2F12.3)") x(i, j), y(i, j)
                enddo
                i = 23
                j = 3
                write (90, "(2I5, 2F12.4)") i, j, x(i, j), y(i, j)
                write (91, "(2F12.3)") x(i, j), y(i, j)
                i = 23
                j = 2
                write (90, "(2I5, 2F12.4)") i, j, x(i, j), y(i, j)
                write (91, "(2F12.3)") x(i, j), y(i, j)
                i = 22
                j = 3
                write (90, "(2I5, 2F12.4)") i, j, x(i, j), y(i, j)
                write (91, "(2F12.3)") x(i, j), y(i, j)
                j = 4
                do i = 21, ib, -1
                    write (90, "(2I5, 2F12.4)") i, j, x(i, j), y(i, j)
                    write (91, "(2F12.3)") x(i, j), y(i, j)
                enddo
                close (90)
                close (91)
                do j = jb, je
                    do i = ib, ie
                        write (12, "(2F12.3)") x(i, j), y(i, j)
                    enddo
                    write (12, "(' A ')")
                enddo
                do i = ib, ie
                    do j = jb, je
                        write (12, "(2F12.3)") x(i, j), y(i, j)
                    enddo
                    write (12, "(' A ')")
                enddo
                j = 3
                do i = 21, 24
                    write (12, "(2F12.3)") x(i, j), y(i, j)
                enddo
                write (12, "(' A ')")
                j = 2
                do i = 22, 23
                    write (12, "(2F12.3)") x(i, j), y(i, j)
                enddo
                write (12, "(' A ')")
                i = 21
                do j = 3, 4
                    write (12, "(2F12.3)") x(i, j), y(i, j)
                enddo
                write (12, "(' A ')")
                i = 22
                do j = 2, 4
                    write (12, "(2F12.3)") x(i, j), y(i, j)
                enddo
                write (12, "(' A ')")
                i = 23
                do j = 2, 4
                    write (12, "(2F12.3)") x(i, j), y(i, j)
                enddo
                write (12, "(' A ')")
                i = 24
                do j = 3, 4
                    write (12, "(2F12.3)") x(i, j), y(i, j)
                enddo
                write (12, "(' A ')")

            endif
            !********** (ntype/=7) endif **********

            do j = jmin, jmax
                do i = imin, imax
                    x(i, j) = x(i, j) - xshift
                    y(i, j) = y(i, j) - yshift
                enddo
            enddo
            if (ntype/=7) then
                write (7, "(1X, 'FINAL X AND Y FIELDS', /)")
                write (7, "(5X, 'I', 10X, 'J', 11X, 'X', 12X, 'Y', /)")
                if (ntype==1) Then
                    do j = jmin, jmax
                        write (9, "(1X, 'J=', I5)") j
                        write (12, "(1X, 'J=', I5)") j
                        do i = imin, imaxo
                            if (ksgi(i,j)==1) Then
                                write (7, "(1X, I10, 5X, I10, 5X, F12.4, 5X, F12.4)") i, j, x(i, j), y(i, j)
                                write (9, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (12, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                    do i = imin, imaxo
                        write (9, "(1X, 'I=', I5)") i
                        write (12, "(1X, 'I=', I5)") i
                        do j = jmin, jmax
                            if (ksgi(i,j)==1) Then
                                write (8, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (12, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                endif
                if (ntype==2) Then
                    do j = jmin, jmaxo
                        write (9, "(1X, 'J=', I5)") j
                        write (12, "(1X, 'J=', I5)") j
                        do i = imin, imax
                            if (ksgi(i,j)==1) Then
                                write (7, "(1X, I10, 5X, I10, 5X, F12.4, 5X, F12.4)") i, j, x(i, j), y(i, j)
                                write (9, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (12, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                    do i = imin, imax
                        write (9, "(1X, 'I=', I5)") i
                        write (12, "(1X, 'I=', I5)") i
                        do j = jmin, jmaxo
                            if (ksgi(i,j)==1) Then
                                write (8, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (12, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                endif
                if (ntype==3) Then
                    do j = jmin, jmax
                        write (9, "(1X, 'J=', I5)") j
                        write (12, "(1X, 'J=', I5)") j
                        do i = imino, imax
                            if (ksgi(i,j)==1) Then
                                write (7, "(1X, I10, 5X, I10, 5X, F12.4, 5X, F12.4)") i, j, x(i, j), y(i, j)
                                write (9, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (12, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                    do i = imino, imax
                        write (9, "(1X, 'I=', I5)") i
                        write (12, "(1X, 'I=', I5)") i
                        do j = jmin, jmax
                            if (ksgi(i,j)==1) Then
                                write (8, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (12, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                endif
                if (ntype==4) Then
                    do j = jmino, jmax
                        write (9, "(1X, 'J=', I5)") j
                        write (12, "(1X, 'J=', I5)") j
                        do i = imin, imax
                            if (ksgi(i,j)==1) Then
                                write (7, "(1X, I10, 5X, I10, 5X, F12.4, 5X, F12.4)") i, j, x(i, j), y(i, j)
                                write (8, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (12, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                    do i = imin, imax
                        write (9, "(1X, 'I=', I5)") i
                        write (12, "(1X, 'I=', I5)") i
                        do j = jmino, jmax
                            if (ksgi(i,j)==1) Then
                                write (8, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (12, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                endif
                if (ntype>=5 .And. ntype<=6) Then
                    do j = jmin, jmax
                        write (9, "(1X, 'J=', I5)") j
                        write (12, "(1X, 'J=', I5)") j
                        do i = imin, imax
                            if (ksgi(i,j)==1) Then
                                write (7, "(1X, I10, 5X, I10, 5X, F12.4, 5X, F12.4)") i, j, x(i, j), y(i, j)
                                write (9, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (12, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                    do i = imin, imax
                        write (9, "(1X, 'I=', I5)") i
                        write (12, "(1X, 'I=', I5)") i
                        do j = jmin, jmax
                            if (ksgi(i,j)==1) Then
                                write (8, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                                write (12, "(1X, F12.4, 5X, F12.4)") x(i, j), y(i, j)
                            endif
                        enddo
                    enddo
                endif
            endif
        else
            open (1, file='gridext.inp', status='unknown')
            do
                read(1, *, IOSTAT=io_status) i, j, x(i, j), y(i, j)
                if (io_status /= 0) then
                    exit
                endif
                x(i, j) = x(i, j) + xshift
                y(i, j) = y(i, j) + yshift
            enddo
            close (1)
        endif
        !****** (ntype/=0) endif ******
        if (isidep==1 .Or. isveg==1) Then
            do j = jmino, jmaxo - 1
                do i = imino, imaxo - 1
                    if (ijct(i,j)>=1 .And. ijct(i,j)<=5) Then
                        dxdi = 0.5*(x(i+1,j+1)-x(i,j+1)+x(i+1,j)-x(i,j))
                        dydi = 0.5*(y(i+1,j+1)-y(i,j+1)+y(i+1,j)-y(i,j))
                        dxdj = 0.5*(x(i+1,j+1)-x(i+1,j)+x(i,j+1)-x(i,j))
                        dydj = 0.5*(y(i+1,j+1)-y(i+1,j)+y(i,j+1)-y(i,j))
                        gii = dxdi*dxdi + dydi*dydi
                        gjj = dxdj*dxdj + dydj*dydj
                        hii = sqrt(gii)
                        hjj = sqrt(gjj)
                        rki(i, j) = hjj/hii
                        rkj(i, j) = hii/hjj
                        xlnutme = 0.25*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1))
                        yltutmn = 0.25*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))
                        hiisq = hii*hii
                        hjjsq = hjj*hjj
                        radsq1 = min(hiisq, hjjsq)
                        radsq2 = max(hiisq, hjjsq)
                        if (isidep==1) Then
                            Call raddep(i, j, xlnutme, yltutmn, radsq1, radsq2, depth)
                            if (depth<=-998.) Then
                                n999 = n999 + 1
                                depcc(i, j) = 0.
                                depfix(i, j) = 0.
                                write (67, "('NO DEP DATA AT I,J,X,Y = ', 2I5, 2E14.5)") i, j, xlnutme, yltutmn
                            else
                                depcc(i, j) = depth
                                depfix(i, j) = depth
                                if (depth>1000.) Then
                                    write (67, "('DEPTH TOO LARGE AT AT I,J,DEPTH = ', 2I5, 2E14.5)") i, j, depth
                                endif
                            endif
                        endif
                        if (isveg==1) Then
                            Call radveg(i, j, xlnutme, yltutmn, radsq1, radsq2, nvtmp)
                            if (nvtmp==0) Then
                                n999 = n999 + 1
                                nvegij(i, j) = 0
                                write (67, "('NO VEG DATA AT I,J,X,Y = ', 2I5, 2E14.5)") i, j, xlnutme, yltutmn
                            else
                                nvegij(i, j) = nvtmp
                            endif
                        endif
                    endif
                enddo
            enddo
        endif
        if (isidep==1) Then
            do j = jmino, jmaxo
                do i = imino, imaxo
                    if (ijct(i,j)>=1 .And. ijct(i,j)<=5) Then
                        rkii(i, j) = 0.5*(rki(i-1,j)+rki(i,j))
                        rkjj(i, j) = 0.5*(rkj(i,j-1)+rkj(i,j))
                        if (ijct(i-1,j)==9) rkii(i, j) = 0.
                        if (ijct(i,j-1)==9) rkjj(i, j) = 0.
                    else
                        rkii(i, j) = 0.
                        rkjj(i, j) = 0.
                    endif
                enddo
            enddo
            do j = jmino, jmaxo
                do i = imino, imaxo
                    xn(i, j) = depcc(i, j)
                enddo
            enddo
            do nd = 1, ndepsm
                do j = jmino, jmaxo
                    do i = imino, imaxo
                        if (ijct(i,j)>=1 .And. ijct(i,j)<=5) Then
                            ctot = rkii(i+1, j) + rkii(i, j) + rkjj(i, j+1) + rkjj(i, j)
                            ccen = (1.-0.0625*ctot)*xn(i, j)
                            cavg = 0.0625*(rkii(i+1,j)*xn(i+1,j)+rkii(i,j)*xn(i-1,j)+rkjj(i,j+1)*xn(i,j+1)+rkjj(i,j)*xn(i,j-1))
                            depcc(i, j) = ccen + cavg
                            if (nd==ndepsm) Then
                                if (ccen<0.) write (67, "('SMOOTH ERR, I,J,CTOT=', 2I5, E14.5)") i, j, ctot
                            endif
                        endif
                    enddo
                enddo
                do j = jmino, jmaxo
                    do i = imino, imaxo
                        if (depfix(i,j)/=0.) depcc(i, j) = depfix(i, j)
                        xn(i, j) = depcc(i, j)
                    enddo
                enddo
            enddo
        endif
        do n = 1, ndepsmf
            do j = jmino, jmaxo
                do i = imino, imaxo
                    xn(i, j) = depcc(i, j)
                enddo
            enddo
            do j = jmino, jmaxo
                do i = imino, imaxo
                    if (ijct(i,j)>=1 .And. ijct(i,j)<=5) Then
                        ctot = rkii(i+1, j) + rkii(i, j) + rkjj(i, j+1) + rkjj(i, j)
                        ccen = (1.-0.0625*ctot)*xn(i, j)
                        cavg = 0.0625*(rkii(i+1,j)*xn(i+1,j)+rkii(i,j)*xn(i-1,j)+rkjj(i,j+1)*xn(i,j+1)+rkjj(i,j)*xn(i,j-1))
                        depcc(i, j) = ccen + cavg
                    endif
                enddo
            enddo
        enddo
        open (39, file='gridext.out', status='unknown')
        do j = jmino, jmaxo
            do i = imino, imaxo
                if (x(i,j)/=0.0 .And. y(i,j)/=0.0) Then
                    write (39, "(2I5, 2X, F10.6, 2X, F10.6)") i, j, x(i, j), y(i, j)
                endif
            enddo
        enddo
        close (39)
        open (30, file='salt.inp', status='unknown')
        write (14, "('C dxdy.inp file, in free format across columns')")
        write (14, "('C')")
        write (14, "('C     I     J        DX            DY           ', 1X, 'DEPTH     BOTTOM ELEV      ZROUGH  VEG TYPE')")
        write (14, "('C')")
        write (16, "('C lxly.inp file, in free format across line')")
        write (16, "('C')")
        write (16, "('C    I     J    XLNUTME       YLTUTMN        CCUE', 1X, '           CCVE          CCUN         CCVN')")
        write (16, "('C')")
        write (15, "(5X, 'I', 5X, 'J', 4X, 'HII', 10X, 'HJJ', 10X, 'HIIHJJ', 6X, 'JACOBIAN', 5X, 'ANG ERROR', /)")
        n999 = 0
        depmax = 0.
        nwcells = 0.
        asqrtg = 0.
        ahihj = 0.
        do j = jmino, jmaxo - 1
            do i = imino, imaxo - 1
                if (ijct(i,j)>=1 .And. ijct(i,j)<=5) Then
                    nwcells = nwcells + 1
                    dxdi = 0.5*(x(i+1,j+1)-x(i,j+1)+x(i+1,j)-x(i,j))
                    dydi = 0.5*(y(i+1,j+1)-y(i,j+1)+y(i+1,j)-y(i,j))
                    dxdj = 0.5*(x(i+1,j+1)-x(i+1,j)+x(i,j+1)-x(i,j))
                    dydj = 0.5*(y(i+1,j+1)-y(i+1,j)+y(i,j+1)-y(i,j))
                    gii = dxdi*dxdi + dydi*dydi
                    gjj = dxdj*dxdj + dydj*dydj
                    gij = dxdi*dxdj + dydi*dydj
                    angd = (acos(gij/sqrt(gii*gjj)))*57.29578
                    angerr = abs(90.-angd)
                    hii = sqrt(gii)
                    hjj = sqrt(gjj)
                    ceu = dxdi/hii
                    cev = dxdj/hjj
                    cnu = dydi/hii
                    cnv = dydj/hjj
                    xlnutme = 0.25*(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1))
                    yltutmn = 0.25*(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))
                    if (isidep==1) Then
                        if (isidptyp==1) Then
                            depth = depcc(i, j)
                            belv = -1.*depth
                            zrough = 0.0
                            depmax = max(depmax, depth)
                        else
                            belv = depcc(i, j)
                            depth = surfelv - belv
                            depth = max(depmin, depth)
                            zrough = 0.0
                            depmax = max(depmax, depth)
                        endif
                    else
                        depth = 0.
                        belv = 0.0
                        zrough = 0.0
                    endif
                    wndshe = 1.0
                    hiihjj = hii*hjj
                    ahihj = ahihj + hiihjj
                    sqrtg = sqrt(gii*gjj-gij*gij)
                    asqrtg = asqrtg + sqrtg
                    hii = hii*hscale
                    hjj = hjj*hscale
                    write (14, "(1X, I5, 2X, I5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, I3)") i, j, hii, hjj, depth, belv, zrough, nvegij(i, j)
                    write (6, "(1X, I5, 2X, I5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, I3)") i, j, hii, hjj, depth, belv, zrough, nvegij(i, j)
                    write (66, "(1X, I5, 2X, I5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, E12.5, 2X, I3)") i, j, hii, hjj, depth, belv, zrough, nvegij(i, j)
                    write (15, "(1X, I5, 1X, I5, 5(2X,E11.4))") i, j, hii, hjj, hiihjj, sqrtg, angerr
                    write (16, "(1X, I5, 1X, I5, 7(1X,E13.6))") i, j, xlnutme, yltutmn, ceu, cev, cnu, cnv, wndshe
                    saltmp = 0.0
                    write (30, "(3I5, 2X, 12F6.1)") nwcells, i, j, saltmp, saltmp, saltmp, saltmp, saltmp, saltmp, saltmp, saltmp, saltmp, saltmp
                    if (ijct(i,j)==1) Then
                        write (27, "(2F12.3, '  -1')") x(i, j), y(i, j)
                        write (27, "(2F12.3, '   1')") x(i+1, j), y(i+1, j)
                        write (27, "(2F12.3, '   1')") x(i, j+1), y(i, j+1)
                        write (27, "(2F12.3, '   1')") x(i, j), y(i, j)
                    endif
                    if (ijct(i,j)==4) Then
                        write (27, "(2F12.3, '  -1')") x(i, j), y(i, j)
                        write (27, "(2F12.3, '   1')") x(i+1, j), y(i+1, j)
                        write (27, "(2F12.3, '   1')") x(i+1, j+1), y(i+1, j+1)
                        write (27, "(2F12.3, '   1')") x(i, j), y(i, j)
                    endif
                    if (ijct(i,j)==3) Then
                        write (27, "(2F12.3, '  -1')") x(i+1, j), y(i+1, j)
                        write (27, "(2F12.3, '   1')") x(i+1, j+1), y(i+1, j+1)
                        write (27, "(2F12.3, '   1')") x(i, j+1), y(i, j+1)
                        write (27, "(2F12.3, '   1')") x(i+1, j), y(i+1, j)
                    endif
                    if (ijct(i,j)==2) Then
                        write (27, "(2F12.3, '  -1')") x(i, j), y(i, j)
                        write (27, "(2F12.3, '   1')") x(i+1, j+1), y(i+1, j+1)
                        write (27, "(2F12.3, '   1')") x(i, j+1), y(i, j+1)
                        write (27, "(2F12.3, '   1')") x(i, j), y(i, j)
                    endif
                    if (ijct(i,j)==5) Then
                        write (27, "(2F12.3, '  -1')") x(i, j), y(i, j)
                        write (27, "(2F12.3, '   1')") x(i+1, j), y(i+1, j)
                        write (27, "(2F12.3, '   1')") x(i+1, j+1), y(i+1, j+1)
                        write (27, "(2F12.3, '   1')") x(i, j+1), y(i, j+1)
                        write (27, "(2F12.3, '   1')") x(i, j), y(i, j)
                    endif
                    if (ijct(i,j)==1) Then
                        write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                        write (25, "('  0', /, 'SEQEND')")
                    endif
                    if (ijct(i,j)==4) Then
                        write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                        write (25, "('  0', /, 'SEQEND')")
                    endif
                    if (ijct(i,j)==3) Then
                        write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                        write (25, "('  0', /, 'SEQEND')")
                    endif
                    if (ijct(i,j)==2) Then
                        write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                        write (25, "('  0', /, 'SEQEND')")
                    endif
                    if (ijct(i,j)==5) Then
                        write (25, "('  0', /, 'POLYLINE', /, '  8', /, 'GRID', /, ' 66', /, '     1')")
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j), y(i+1, j)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i+1, j+1), y(i+1, j+1)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j+1), y(i, j+1)
                        write (25, "('  0', /, 'VERTEX', /, '  8', /, 'GRID', /, ' 10', /, F12.3, /, ' 20', /, F12.3)") x(i, j), y(i, j)
                        write (25, "('  0', /, 'SEQEND')")
                    endif
                    icomp(nwcells) = i
                    jcomp(nwcells) = j
                    xcomp(nwcells) = xlnutme
                    ycomp(nwcells) = yltutmn
                endif
            enddo
        enddo
        aerr = 2.*abs(asqrtg-ahihj)/(asqrtg+ahihj)
        write (15, "(//)")
        write (15, "(1X, 'ASQRTG=', 2X, E11.4, 2X, 'ASHIHJ=', 2X, E11.4, 2X, 'AERR=', 2X, E11.4)") asqrtg, ahihj, aerr
        write (15, "(1X, 'NWCELLS=', I10)") nwcells
        write (6, "(1X, 'NWCELLS=', I10)") nwcells
        write (66, "(1X, 'NWCELLS=', I10)") nwcells
        write (6, "(1X, 'N999 =', I10)") n999
        write (66, "(1X, 'N999 =', I10)") n999
        close (30)
        if (isgg/=0) then
            open (17, file='gcell.inp', status='unknown')
            open (18, file='gcellmap.out', status='unknown')
            write (18, "('C gcellmap.inp file, in free format across columns')")
            write (18, "('C')")
            write (18, "('  IGRAPHIC  JGRAPHIC     ICOMP     JCOMP         RMIN')")
            write (18, "('C')")
            do is = 1, 4
                read (17, "(130X)")
            enddo
            read (17, "(I3, 2X, 125I1)") jctmp
            close (17)
            open (17, file='gcell.inp', status='unknown')
            do is = 1, 4
                read (17, "(130X)")
            enddo
            if (jctmp/=jgm) Then
                do jt = 1, jgm, 120
                    jf = jt
                    jlast = jt + 119
                    if (jlast>jgm) jlast = jgm
                    write (7, "(1X, '  GCELL TYPE ARRAY,J=', I5, 2X, 'TO J=', I5, //)") jf, jlast
                    do i = 1, igm
                        read (17, "(120I1)")(ijctg(i,j), j=jf, jlast)
                        write (7, "(1X, I3, 2X, 125I1)") i, (ijctg(i,j), j=jf, jlast)
                    enddo
                enddo
            else
                do it = 1, igm, 120
                    ifirst = it
                    ilast = it + 119
                    if (ilast>igm) ilast = igm
                    write (7, "(1X, '  GELL TYPE ARRAY,I=', I5, 2X, 'TO I=', I5, //)") ifirst, ilast
                    do j = jgm, 1, -1
                        read (17, "(I3, 2X, 125I1)") jdumy, (ijctg(i,j), i=ifirst, ilast)
                        write (7, "(1X, I3, 2X, 125I1)") jdumy, (ijctg(i,j), i=ifirst, ilast)
                    enddo
                enddo
            endif
            nwgg = 0
            write (18, "(4I10, 5X, F12.4)") igm, jgm, nwtgg
            do jg = 1, jgm
                write (6, "(' JG = ', I5)") jg
                write (66, "(' JG = ', I5)") jg
                do ig = 1, igm
                    if (ijctg(ig,jg)/=0) Then
                        nwgg = nwgg + 1
                        write (6, "(' NWGG = ', I5)") nwgg
                        write (66, "(' NWGG = ', I5)") nwgg
                        xgg = cdlon1 + (cdlon2*float(ig)+cdlon3)/60.
                        ygg = cdlat1 + (cdlat2*float(jg)+cdlat3)/60.
                        if (nwtgg>=1) Then
                            rmin1 = 1.E+9
                            do nw = 1, nwcells
                                rtmp = sqrt((xcomp(nw)-xgg)**2+(ycomp(nw)-ygg)**2)
                                if (rtmp<rmin1) Then
                                    nw1 = nw
                                    rmin1 = rtmp
                                    icompt1 = icomp(nw)
                                    jcompt1 = jcomp(nw)
                                endif
                            enddo
                        endif
                        if (nwtgg>=2) Then
                            rmin2 = 1.E+9
                            do nw = 1, nwcells
                                if (nw/=nw1) Then
                                    rtmp = sqrt((xcomp(nw)-xgg)**2+(ycomp(nw)-ygg)**2)
                                    if (rtmp<rmin2) Then
                                        nw2 = nw
                                        rmin2 = rtmp
                                        icompt2 = icomp(nw)
                                        jcompt2 = jcomp(nw)
                                    endif
                                endif
                            enddo
                        endif
                        if (nwtgg>=3) Then
                            rmin3 = 1.E+9
                            do nw = 1, nwcells
                                if (nw/=nw1 .And. nw/=nw2) Then
                                    rtmp = sqrt((xcomp(nw)-xgg)**2+(ycomp(nw)-ygg)**2)
                                    if (rtmp<rmin3) Then
                                        nw3 = nw
                                        rmin3 = rtmp
                                        icompt3 = icomp(nw)
                                        jcompt3 = jcomp(nw)
                                    endif
                                endif
                            enddo
                        endif
                        if (nwtgg>=4) Then
                            rmin4 = 1.E+9
                            do nw = 1, nwcells
                                if (nw/=nw1 .And. nw/=nw2 .And. nw/=nw3) Then
                                    rtmp = sqrt((xcomp(nw)-xgg)**2+(ycomp(nw)-ygg)**2)
                                    if (rtmp<rmin4) Then
                                        nw4 = nw
                                        rmin4 = rtmp
                                        icompt4 = icomp(nw)
                                        jcompt4 = jcomp(nw)
                                    endif
                                endif
                            enddo
                        endif
                        if (nwtgg>=1) write (18, "(4I10, 5X, F12.4)") ig, jg, icompt1, jcompt1, rmin1
                        if (nwtgg>=2) write (18, "(4I10, 5X, F12.4)") ig, jg, icompt2, jcompt2, rmin2
                        if (nwtgg>=3) write (18, "(4I10, 5X, F12.4)") ig, jg, icompt3, jcompt3, rmin3
                        if (nwtgg>=4) write (18, "(4I10, 5X, F12.4)") ig, jg, icompt4, jcompt4, rmin4
                    endif
                enddo
            enddo
            write (18, "(4I10, 5X, F12.4)") nwgg
            close (17)
            close (18)
        endif
    endif
    !********** (ntype>=8) endif **********
    write (6, "(' DEPMAX = ', E12.5)") depmax
    write (66, "(' DEPMAX = ', E12.5)") depmax
    close (2)
    close (3)
    close (4)
    close (7)
    close (8)
    close (9)
    close (10)
    close (12)
    close (13)
    close (14)
    close (15)
    close (16)
    close (66)
    close (67)
    write (25, "('  8', /, 'GRID', /, '  0', /, 'ENDSEC', /, '  0', /, 'EOF')")
    write (26, "('  8', /, 'GRID', /, '  0', /, 'ENDSEC', /, '  0', /, 'EOF')")
    close (25)
    close (26)
    close (27)
    close (89)
    deallocate(depd, depe, depn, nvegd, vege, vegn, nnveg, x, y)

    end Program gefdc


    subroutine radveg(i, j, utme, utmn, radsq1, radsq2, nvtmp)
    use OneModule
    use TwoModule
    use VegModule
    real(kind=8)::utme, utmn, radsq1, radsq2, nvtmp, radsq
    integer::i, j, ndum, ntmp
    parameter (nddmax=225000, nvdmax=1000)

    radsq = radm*radm*radsq1
    do ntmp = 0, nvegtyp
        nnveg(ntmp) = 0
    enddo
    do n = 1, nvegdat
        rstmp = (vege(n)-utme)**2 + (vegn(n)-utmn)**2
        if (rstmp<=radsq) Then
            nttmp = nvegd(n)
            nnveg(nttmp) = nnveg(nttmp) + 1
        endif
    enddo
    write (6, '(16I5)') i, j, (nnveg(n), n=1, nvegtyp)
    write (66, '(16I5)') i, j, (nnveg(n), n=1, nvegtyp)
    ndum = -1
    do ntmp = 1, nvegtyp
        if (nnveg(ntmp)>ndum) Then
            ndum = nnveg(ntmp)
            nvtmp = ntmp
        endif
    enddo
    if (nnveg(nvtmp)==0) nvtmp = 0
    return
    end subroutine radveg



    subroutine raddep(i, j, utme, utmn, radsq1, radsq2, depth)
    use OneModule
    use TwoModule
    use CxyModule
    real(kind=8):: sumw, sumwd, radsq, wt, adep, xsw, ysw, xse, yse, xne, yne, xnw, ynw, a11, a12, a13, a21, a22, a23, a31, a32, a33,  detb, cixx, ciyy, cixy, cjxx, cjyy, cjxy
    integer:: ndum, ntmp, ip ,jp
    real(kind=8) :: utme, utmn, radsq1, radsq2
    integer :: i, j
    if (radm<0.) then
        radsq1 = radm*radm*radsq1
        radsq2 = radm*radm*radsq2
        radsq = max(radsq1, radsq2)
        xsw = 0.
        ysw = 0.
        xse = x(i+1, j) - x(i, j)
        yse = y(i+1, j) - y(i, j)
        xne = x(i+1, j+1) - x(i, j)
        yne = y(i+1, j+1) - y(i, j)
        xnw = x(i, j+1) - x(i, j)
        ynw = y(i, j+1) - y(i, j)
        a11 = xse
        a12 = yse
        a13 = a11*a12
        a21 = xne
        a22 = yne
        a23 = a21*a22
        a31 = xnw
        a32 = ynw
        a33 = a31*a32
        detb = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a31*a22*a13 - a32*a23*a11 - a33*a21*a12
        detb = 1./detb
        cixx = (a22*a33+a13*a32-a32*a23-a33*a12)*detb
        ciyy = (a11*a33+a23*a31-a31*a13-a33*a21)*detb
        cixy = (a12*a31+a21*a32-a31*a22-a32*a11)*detb
        cjxx = (a12*a23+a13*a32-a22*a13-a33*a12)*detb
        cjyy = (a11*a33+a13*a21-a31*a13-a23*a11)*detb
        cjxy = (a11*a22+a12*a31-a32*a11-a21*a12)*detb
        sumw = 0.
        sumwd = 0.
        ip = i + 1
        jp = j + 1
        depmxx = -1000.
        depmnn = 1000.
        do n = 1, ndepdat
            xval = depe(n) - x(i, j)
            yval = depn(n) - y(i, j)
            ival = i + cixx*xval + ciyy*yval + cixy*xval*yval
            jval = j + cjxx*xval + cjyy*yval + cjxy*xval*yval
            if (ival>=i .And. ival<=ip) Then
                if (jval>=j .And. jval<=jp) Then
                    if (cdep>0.001) Then
                        rstmp = (depe(n)-utme)**2 + (depn(n)-utmn)**2
                        if (rstmp<radsq) Then
                            if (i>=12 .And. j<=40) Then
                                write (89, "(2I5, 4F12.2, F10.2)") i, j, utme, utmn, depe(n), depn(n), depd(n)
                            endif
                            wt = (radsq-rstmp)/(radsq+rstmp)
                            wt = wt**cdep
                            sumw = sumw + wt
                            sumwd = sumwd + wt*depd(n)
                        endif
                    else
                        if (i>=12 .And. j<=40) Then
                            write (89, "(2I5, 4F12.2, F10.2)") i, j, utme, utmn, depe(n), depn(n), depd(n)
                        endif
                        wt = 1.
                        sumw = sumw + wt
                        sumwd = sumwd + wt*depd(n)
                    endif
                    depmxx = max(depmxx, depd(n))
                    depmnn = min(depmnn, depd(n))
                endif
            endif
        enddo
        if (sumw==0.) Then
            depth = -999.
        else
            depth = sumwd/sumw
            if (isidptyp==2 .And. depth>surfelv) depth = -999.
            write (89, "('SUMMARY I,J,DMN,DAV,DMX = '2I5, 3F10.2)") i, j, depmnn, depth, depmxx
            adep = abs(depth)
            if (adep>1000.) Then
                write (67, "('ERR AT I,J,SUMW,SUMWD,RADSQ,RSTMP=', 2I5, 4E14.5)") i, j, sumw, sumwd, radsq, rstmp
            endif
        endif
    else
        sumw = 0.
        sumwd = 0.
        radsq = radm*radm*radsq1
        do n = 1, ndepdat
            rstmp = (depe(n)-utme)**2 + (depn(n)-utmn)**2
            if (rstmp<radsq) Then
                wt = (radsq-rstmp)/(radsq+rstmp)
                wt = wt**cdep
                sumw = sumw + wt
                sumwd = sumwd + wt*depd(n)
            endif
        enddo
        if (sumw==0.) Then
            depth = -999.
        else
            depth = sumwd/sumw
            adep = abs(depth)
            if (adep>1000.) Then
                write (67, "('ERR AT I,J,SUMW,SUMWD,RADSQ,RSTMP=', 2I5, 4E14.5)") i, j, sumw, sumwd, radsq, rstmp
            endif
        endif
    endif
    return
    end subroutine raddep



    subroutine vutmbay
    use CutmModule
    real(kind=8):: er, pr
    rk = .9996
    phi0 = 0.0
    cmerid = 75.0
    er = 6378206.4D0
    pr = 6356583.8D0
    ae = er
    ecc2 = 1. - pr**2./er**2
    ecc2p = ecc2/(1.-ecc2)
    return
    end subroutine vutmbay


    subroutine calsmd
    use CsmdModule
    use CxyModule
    real(kind=8) :: sxtmp, ci, cip1, cip2, cim1, cim2, cj, cjp1, cjp2
    if (ijsmd/=0) Then
        sxtmp = 0.5*(x(ie,jb)-x(ib,jb)+x(ie,je)-x(ib,je))
        do j = jb + 1, je - 1
            sxtmp = sxtmp + x(ie, j) - x(ib, j)
        enddo
        sytmp = 0.5*(y(ib,je)-y(ib,jb)+y(ie,je)-y(ie,jb))
        do i = ib + 1, ie - 1
            sytmp = sytmp + y(i, je) - y(i, jb)
        enddo
        smd = sxtmp/sytmp
        return
    endif
    if (ismd/=0) Then
        if (ismd==ib) Then
            ci = 0.5*(x(ismd,jb)+x(ismd,je))
            cip1 = 0.5*(x(ismd+1,jb)+x(ismd+1,je))
            cip2 = 0.5*(x(ismd+2,jb)+x(ismd+2,je))
            do j = jb + 1, je - 1
                ci = ci + x(ismd, j)
                cip1 = cip1 + x(ismd+1, j)
                cip2 = cip2 + x(ismd+2, j)
            enddo
            sxtmp = -1.5*ci + 2.*cip1 - 0.5*cip2
        endif
        if (ismd>ib .And. ismd<ie) Then
            cim1 = 0.5*(x(ismd-1,jb)+x(ismd-1,je))
            cip1 = 0.5*(x(ismd+1,jb)+x(ismd+1,je))
            do j = jb + 1, je - 1
                cim1 = cim1 + x(ismd-1, j)
                cip1 = cip1 + x(ismd+1, j)
            enddo
            sxtmp = 0.5*(cip1-cim1)
        endif
        if (ismd==ie) Then
            ci = 0.5*(x(ismd,jb)+x(ismd,je))
            cim1 = 0.5*(x(ismd-1,jb)+x(ismd-1,je))
            cim2 = 0.5*(x(ismd-2,jb)+x(ismd-2,je))
            do j = jb + 1, je - 1
                ci = ci + x(ismd, j)
                cim1 = cim1 + x(ismd-1, j)
                cim2 = cim2 + x(ismd-2, j)
            enddo
            sxtmp = 1.5*ci - 2.*cim1 + 0.5*cim2
        endif
        sytmp = y(ismd, je) - y(ismd, jb)
        smd = sxtmp/sytmp
        return
    endif
    if (jsmd/=0) Then
        sxtmp = x(ie, jsmd) - x(ib, jsmd)
        if (jsmd==jb) Then
            cj = 0.5*(y(ib,jsmd)+y(ie,jsmd))
            cjp1 = 0.5*(y(ib,jsmd+1)+y(ie,jsmd+1))
            cjp2 = 0.5*(y(ib,jsmd+2)+y(ie,jsmd+2))
            do i = ib + 1, ie - 1
                cj = cj + y(i, jsmd)
                cjp1 = cjp1 + y(i, jsmd+1)
                cjp2 = cjp2 + y(i, jsmd+2)
            enddo
            sytmp = -1.5*cj + 2.*cjp1 - 0.5*cjp2
        endif
        if (jsmd>jb .And. jsmd<je) Then
            cjm1 = 0.5*(y(ib,jsmd-1)+y(ie,jsmd-1))
            cjp1 = 0.5*(y(ib,jsmd+1)+y(ie,jsmd+1))
            do i = ib + 1, ie - 1
                cjm1 = cjm1 + y(i, jsmd-1)
                cjp1 = cjp1 + y(i, jsmd+1)
            enddo
            sytmp = 0.5*(cjp1-cjm1)
        endif
        if (jsmd==je) Then
            cj = 0.5*(y(ib,jsmd)+y(ie,jsmd))
            cjm1 = 0.5*(y(ib,jsmd-1)+y(ie,jsmd-1))
            cjm2 = 0.5*(y(ib,jsmd-2)+y(ie,jsmd-2))
            do i = ib + 1, ie - 1
                cj = cj + y(i, jsmd)
                cjm1 = cjm1 + y(i, jsmd-1)
                cjm2 = cjm2 + y(i, jsmd-2)
            enddo
            sytmp = 1.5*cj - 2.*cjm1 + 0.5*cjm2
        endif
        smd = sxtmp/sytmp
        return
    endif
    return
    end subroutine calsmd

