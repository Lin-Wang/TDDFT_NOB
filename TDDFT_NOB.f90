            subroutine NOB_select()
  ! This code handle the collapsing of natural orbitals
            implicit none
            integer iislda,kpt
            integer i1,i2,ierr
            real*8 xx,yy,EE0,EEnew
            integer, allocatable, dimension(:) :: list_tmp,list_tmp2
            real*8, external :: ran1
            integer indx, ind_try,num_select
            real*8 pi,dn,sum1
            complex*16,allocatable,dimension(:) :: zdn1,zdn2
            real*8, allocatable,dimension(:) :: prob1,prob2,aprob1,aprob2
            integer, allocatable,dimension(:) :: iselect0_tmp,iselect_tmp
            integer, allocatable,dimension(:) :: iflag,iflag2,ind_tmp
            real*8 delta,sum2
            integer mm0,nn0,num0,num2,mm,nn,num,ntry_tmp


            allocate(list_tmp(mst_win))
            allocate(list_tmp2(mst_win))

            do i1=1,20
             xx=ran1(irandom_tddft)
             enddo


 2022      continue

           num_try=num_try+1
         
            goto 2021
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc This is the old selection, which is not quite correct
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            do iislda=1,islda
            do kpt=1,nkpt

            list_tmp=0
            list_tmp2=0
            num_select=0
 3022          continue               
               xx=ran1(irandom_tddft)
               ind_try=xx*mst_win+1
               if(ind_try.gt.mst_win) ind_try=mst_win


               if(list_tmp(ind_try).eq.0) then ! This state is not taken yet
                 xx=ran1(irandom_tddft)
                 if(xx.lt.alambda_td(ind_try,kpt,iislda)) then  ! select this one
                 list_tmp(ind_try)=1
                 num_select=num_select+1
                 list_tmp2(num_select)=ind_try
                 if(num_select.ge.numb_td(kpt,iislda)) goto 3023   ! jump out
                 endif
               endif
               goto 3022
 3023           continue
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
             list_select(1:num_select,kpt,iislda)=list_tmp2(1:num_select)
             enddo  ! kpt
             enddo  ! iislda
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc This is the old selection, which is not quite correct
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccc JUMP out point for the old version
2021         continue
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc The following is the new selection method
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            allocate(iselect0_tmp(mst_win))
            allocate(iselect_tmp(mst_win))
            allocate(iflag(mst_win))
            allocate(iflag2(mst_win))
            allocate(ind_tmp(mst_win))
            allocate(prob1(mst_win))
            allocate(prob2(mst_win))
            allocate(aprob1(mst_win))
            allocate(aprob2(mst_win))


            do iislda=1,islda
            do kpt=1,nkpt

            mm0=mst_win
            nn0=numb_td(kpt,iislda)

            delta=0.01
            if(delta.gt.1.d0/(nn0+1)) delta=1.d0/(nn0+1) 
            ! So, not matter what, I cannot select num0 larger than nn0

2003        continue
            num0=0
            num2=0
            iselect0_tmp=0
            iflag=0
            iflag2=0
            do i=1,mm0
            if(alambda_td(i,kpt,iislda).gt.1.d0-delta) then
            num2=num2+1
            iflag2(i)=1
            xx=ran1(irandom_tddft)
            if(xx.lt.alambda_td(i,kpt,iislda)) then
            num0=num0+1
            iselect0_tmp(num0)=i
            iflag(i)=1
            endif
            endif
            enddo

2004        continue
!ccccccccc num2 is preconsidered, so not included in the following SUMFORD
!procedure, and num0 is actually selected, so nn0 should be reduced accordingly

            mm=mm0-num2
            nn=nn0-num0
            if(nn.gt.mm) goto 2003   ! unselected too many large lambda

            num=0
            do i=1,mm0
            if(iflag2(i).eq.0) then   ! considered before
            num=num+1
            ind_tmp(num)=i
            prob1(num)=alambda_td(i,kpt,iislda)
            endif
            enddo

            if(num.ne.mm) then
            write(6,*) "SOMETHING WRONG in NOB Selection, stop"
            stop
            endif

            sum1=0.d0
            do i=1,mm
            sum1=sum1+prob1(i)
            enddo

            do i=1,mm
            prob1(i)=prob1(i)/sum1
            if(nn*prob1(i).ge.1.d0-1.E-5) then ! also preselect this
            num0=num0+1    ! num0 is the number actually preselected
            num2=num2+1    ! num2 is the number, considered in preelection, so  not in the Sumford procedure
            iselect0_tmp(num0)=ind_tmp(i)
            iflag(ind_tmp(i))=1     ! the flag this i is actually preselected
            iflag2(ind_tmp(i))=1    ! the flag this i is considered in preselection
            goto 2004   ! preselect this i, not in the main Sumford selection procedure
            endif
            enddo
!ccccccccccccccccc  passed all the preselect here

            sum2=0.d0
            do i=1,mm
            prob2(i)=prob1(i)/abs(1-nn*prob1(i))
            sum2=sum2+prob2(i)
            enddo
            do i=1,mm
            prob2(i)=prob2(i)/sum2
            enddo

            do i=1,mm
            if(i.eq.1) then
            aprob1(i)=prob1(i)
            aprob2(i)=prob2(i)
            else
            aprob1(i)=aprob1(i-1)+prob1(i)
            aprob2(i)=aprob2(i-1)+prob2(i)
            endif
            enddo
            aprob1(mm)=1.00001d0
            aprob2(mm)=1.00001d0

            ntry_tmp=0
2005        continue
            ntry_tmp=ntry_tmp+1

            iflag(1:mm)=0
            num=0
            xx=ran1(irandom_tddft)
            do i=1,mm
            if(xx.lt.aprob1(i)) then  ! this is the first round selection
            num=num+1
            iselect_tmp(num)=i
            iflag(i)=1
            goto 2006
            endif
            enddo
2006        continue

            do iii=1,nn-1
            xx=ran1(irandom_tddft)
            do i=1,mm
            if(xx.lt.aprob2(i)) then
              if(iflag(i).eq.1) goto 2005  ! double select, abandon it
            num=num+1
            iselect_tmp(num)=i
            iflag(i)=1
            goto 2007
            endif
            enddo
2007        continue
            enddo
!  Now, here we have successfully choosed one nn lambda
            if(inode_tot.eq.1) then
            write(6,"('ntry(SUMFORD),nn,kpt,islda', 3(i4,1x))")  ntry_tmp,nn,kpt,iislda
            endif
             do i=1,num0 
             list_select(i,kpt,iislda)=iselect0_tmp(i)
             enddo
             do i=1,nn
             list_select(num0+i,kpt,iislda)=ind_tmp(iselect_tmp(i))
             enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             enddo   ! kpt
             enddo   ! iislda
       
            deallocate(iselect0_tmp)
            deallocate(iselect_tmp)
            deallocate(iflag)
            deallocate(iflag2)
            deallocate(ind_tmp)
            deallocate(prob1)
            deallocate(prob2)
            deallocate(aprob1)
            deallocate(aprob2)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc This is  the new selection method
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       


             EE0=0.d0
             EEnew=0.d0
             do iislda=1,islda
             do kpt=1,nkpt

             do i1=1,mst_win
             EE0=EE0+E_td(i1,kpt,iislda)*alambda_td(i1,kpt,iislda)* &
                 2.d0/islda*weighkpt_2(kpt)
             enddo

             do i2=1,numb_td(kpt,iislda)
             i1=list_select(i2,kpt,iislda)
             EEnew=EEnew+E_td(i1,kpt,iislda)*2.d0/islda*weighkpt_2(kpt)
             enddo

             enddo
             enddo
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccc

             if(iselect_nob_opt.eq.2.or.iselect_nob_opt.eq.3) goto 2023  ! always pass from the first choice here
          
             if(iselect_nob_opt.eq.1) then
             if(EEnew.lt.EE0) then
            ! successful, return
             goto 2023
             else
             yy=exp(-(EEnew-EE0)/beta_kT_tddft) 
             xx=ran1(irandom_tddft)
             if(xx.lt.yy) then
             goto 2023
             endif
             endif
             endif

             if(iselect_nob_opt.ne.1.and. &
                iselect_nob_opt.ne.2.and. &
                iselect_nob_opt.ne.3) then
             write(6,*) "iselect_nob_opt must be: 1,2,3 for TDDFT_NOB,stop"
             stop
             endif
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
           ! Not successful, go back
             goto 2022

 2023        continue

             if(inode_tot.eq.1.and.iselect_nob_opt.eq.1) then
             write(6,"('NOB select, num_try,EE0,EEnew=',i4,2(E15.7,1x))"), num_try,EE0*27.211396,EEnew*27.211396
             endif
            deallocate(list_tmp)
            deallocate(list_tmp2)


            call mpi_bcast(list_select,mst_win*nkpt*islda,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            ! To make sure all the nodes have chosen the same eigen states

            select_weight=0.d0
            do iislda=1,islda
            do kpt=1,nkpt
            do i1=1,numb_td(kpt,iislda)
            i2=list_select(i1,kpt,iislda)
            select_weight(i2,kpt,iislda)=1.d0
            enddo
            enddo
            enddo

            !cccccccccccccccccccccccccccccccccccccc

            pi=4*datan(1.d0)
            allocate(zdn1(mst_win))
            allocate(zdn2(mst_win))

            do iislda=1,islda
            do kpt=1,nkpt
            zdn1=cmplx(0.d0,0.d0)
            zdn2=cmplx(0.d0,0.d0)
            do i1=1,mst_win
            dn=select_weight(i1,kpt,iislda)-alambda_td(i1,kpt,iislda)
            xx=ran1(irandom_tddft)
            if(dn.gt.0.d0) then
              zdn1(i1)=abs(dn)*exp(-2*pi*xx*cmplx(0.d0,1.d0))   ! random phase
            else
              zdn2(i1)=abs(dn)*exp(-2*pi*xx*cmplx(0.d0,1.d0))
            endif
            enddo

            TCD_alambda(:,:,kpt,iislda)=cmplx(0.d0,0.d0)
            do i1=1,mst_win
            do i2=1,mst_win
            if(i1.ne.i2) then
            TCD_alambda(i1,i2,kpt,iislda)=conjg(zdn1(i1))*zdn2(i2)+zdn1(i2)*conjg(zdn2(i1))
            endif
            enddo
            enddo


  !  TCD_almabda is the matrix based on cpsi_td (Natural orbital), 
  !  it is used to generate TCD_All, then to calculate fatom_TCD
  !  fatom_TCD can be used to scale the velocity in that direction (1 degree of
  !  freedom). It can also be used to decide whether to take this collapsing
  !  (whether this degree of freedom has enough kinetic energy). 

  ! Since the random number is a bit dangerous, so bcast again


            if(inode_tot.eq.1) then
              write(6,*) "INSIDE SELECT, TRY TO COLLAPSE"
              write(6,*) "Try: alambda collapse: kpt,iislda,num_try",kpt,iislda,num_try
            do i1=1,mst_win
            write(6,"('i1,almabda(old),alambda(new) ',i4,2(f17.14,1x),1x,E14.7)")  &
            i1,alambda_td(i1,kpt,iislda),select_weight(i1,kpt,iislda),E_td(i1,kpt,iislda)
            enddo
            endif

            enddo
            enddo
            deallocate(zdn1)
            deallocate(zdn2)


            call mpi_bcast(TCD_alambda,mst_win*mst_win*nkpt*islda,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

            
             return
            end  subroutine NOB_select



            subroutine NOB_proj3()
    ! This code prepare for iproj_occ=3 style calculation
            implicit none
            integer iislda,kpt,m,m1

            do iislda=1,islda
            do kpt=1,nkpt
            call ugIOBP(cpsi_td(1,1,kpt,iislda),kpt,1,5,iislda,-1,nkpt,islda)
            ! store the wave function in 5, to be used for iproj_occ=3 stype
            ! Etotcalc.f calculation
            enddo
            enddo

            do iislda=1,islda
            do kpt=1,nkpt
            do m=1,mx 
            imap_proj_sp(m,kpt,iislda)=m   ! probably not used for this
            imap_proj(m,kpt,iislda)=m
            enddo

            occ_proj_sp(:,kpt,iislda)=0.d0
            do m=1,mst_win0-1
            occ_proj_sp(m,kpt,iislda)=weighkpt_2(kpt)*2.d0/islda
            enddo

            do m=1,numb_td(kpt,iislda)
            m1=list_select(m,kpt,iislda)+mst_win0-1
            occ_proj_sp(m1,kpt,iislda)=weighkpt_2(kpt)*2.d0/islda
            enddo
            enddo
            enddo


            do iislda=1,islda
            if(iislda.eq.1) then
            iferup_iproj3=0.d0 
            do kpt=1,nkpt
            do m=1,mst_win0-1
            iferup_iproj3(m,kpt)=1   ! this is a special flag, to deal with the special 1 to 1 mapping
            enddo
            enddo
            elseif(iislda.eq.2) then
            iferdw_iproj3=0.d0 
            do kpt=1,nkpt
            do m=1,mst_win0-1
            iferdw_iproj3(m,kpt)=1
            enddo
            enddo
            endif
            enddo

            return
            end subroutine NOB_proj3




            subroutine NOB_cdens()
       ! This code prepare for cdens from the coeff_proj_occ in Etot
       ! calculation
            implicit none
            integer iislda,kpt,i1,i2
            complex*16, allocatable, dimension(:,:) :: coeff_tmp
            complex*16 one,zero
            complex*16,allocatable,dimension(:,:) :: dens_tmp
            complex*16,allocatable,dimension(:) :: workx
            real*8,allocatable,dimension(:) :: workrx
            integer lwork,info
            real*8 sum1

            allocate(coeff_tmp(mst_win,mst_win))

            do iislda=1,islda
            do kpt=1,nkpt

            do i1=1,mst_win
            do i2=1,mst_win
            coeff_tmp(i1,i2)=coeff_proj_occ(i1+mst_win0-1,i2+mst_win0-1,kpt,iislda)*select_weight(i2,kpt,iislda)
            enddo
            enddo


            if(inode_tot.eq.1) then
            write(6,*) "alambda collapse: kpt,iislda,num_try",kpt,iislda,num_try
            do i1=1,mst_win
            write(6,"('i1,almabda(old),alambda(new) ',i4,2(f17.14,1x),1x,E14.7)")  &
            i1,alambda_td(i1,kpt,iislda),select_weight(i1,kpt,iislda),E_td(i1,kpt,iislda)
            enddo
            endif
            do i1=1,mst_win
            alambda_td(i1,kpt,iislda)=select_weight(i1,kpt,iislda)
            enddo


            one=cmplx(1.d0,0.d0)
            zero=cmplx(0.d0,0.d0)

            call zgemm('n','c',mst_win,mst_win,mst_win,one,coeff_proj_occ(mst_win0,mst_win0,kpt,iislda),mst,coeff_tmp,mst_win, &
             zero,cdens_td(1,1,kpt,iislda),mst_win)
    ! should update cdens_td, instead of cdens_td, since cdens_td0=cdesen_td
    ! will be assigned after this step, in TDDFT_OUTPUT.f/bak_tddft.f call
    !ccccccccccccccccccccccccccccccccccccccccccc
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if(inode_tot.eq.1) then
            write(6,*) "Dens(i1,i1),Dens_offdiag in eigenstate basis"
            do i1=1,mst_win
            sum1=0.d0
            do i2=1,mst_win
            if(i2.ne.i1) then
            sum1=sum1+abs(cdens_td(i2,i1,kpt,iislda))+abs(cdens_td(i1,i2,kpt,iislda))
            endif
            enddo
            write(6,"('i1, Dens_diag,Dens_off ',i4,2(E14.7,1x))")  &
            i1,abs(cdens_td(i1,i1,kpt,iislda)),sum1
            enddo
            endif


            enddo
            enddo
            deallocate(coeff_tmp)
            return
            end subroutine NOB_cdens
