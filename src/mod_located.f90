module located
        contains      
                subroutine locate(iunit,string)
                integer :: ii        
                integer iunit
                character string*(*)
                character*80 linia

                rewind(iunit)

                ii=0
               do while(ii.eq.0)
                 read(iunit,"(a80)", end=10)linia
                 if(index(linia,string).ne.0) then
                 ii=1
                 return
                !   found=.true.
                 end if
               end do
                10 write(*,*) string, 'section not found'
                return
                end
end module located

