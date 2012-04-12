##################################################################################################
#
# 	checkdeallocations.sh searches all source files for "allocate(...)" statements. The output
# 	from grep is trimmed to reveal the variable to which memory is allocated in the code, so that
# 	a corresponding deallocate statement may be searched for. 
#
#	Please note that this is most certainly NOT a thorough check. Deallocations are called in 
#	many different subroutines, which this script does NOT consider. Furthermore, it doesn't
#	check for deallocations that have been "commented out".
#
##################################################################################################
echo ''
echo 'The following arrays are not deallocated explicitly:'
echo '----------------------------------------------------'
for i in `grep -oh '\ballocate(\b'.*\) *.f90 | sed 's/allocate(//' | sed 's/(.*//' | sed 's/)//'`
do
	grep -q deallocate\($i\) *.f90 && : || echo '    '$i 
done
echo ''
