program identy
	implicit none
	real(8),dimension(5,5)  :: Id
	integer(4)		:: i,j

Id = 0.0d0
do j = 1, 4
	do i = 1, 4
	if(i == j) Id(i,j) = 1.0d0 
	end do
end do

do j = 1, 4
	do i = 1, 4
		print *, i,j,Id(i,j)
	end do
end do
end program identy
