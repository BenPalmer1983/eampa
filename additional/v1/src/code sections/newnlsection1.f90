!----------------------------------------------------  
! Build Inner Volume
!---------------------------------------------------- 	  
	  m = 0
!loop over copies to make full crystal	  
	  coordsStart = configHeaderI(i,headerWidth-1)
	  coordsLength = configHeaderI(i,headerWidth)	  
	  do x=1,xCopy
	    do y=1,yCopy
	      do z=1,zCopy
		    do n=coordsStart,(coordsStart+coordsLength-1)  
			  m = m + 1
			  atomCounter = atomCounter + 1
!get co-ords with lattice parameter and copy applied
			  innerKey(m) = configCoordsI(n)
			  xCoord = alat * (x + configCoordsR(n,1) - 1)
			  yCoord = alat * (y + configCoordsR(n,2) - 1)
			  zCoord = alat * (z + configCoordsR(n,3) - 1)
!Apply configuration unit vector to these co-ordinates
!x-coord
              innerCoords(m,1) = xCoord*workingUnitVector(1,1) + &
              yCoord*workingUnitVector(1,2) + zCoord*workingUnitVector(1,3)
!y-coord
              innerCoords(m,2) = xCoord * workingUnitVector(2,1) + &
              yCoord * workingUnitVector(2,2) + zCoord * workingUnitVector(2,3)
!z-coord
              innerCoords(m,3) = xCoord * workingUnitVector(3,1) + &
              yCoord * workingUnitVector(3,2) + zCoord * workingUnitVector(3,3)
!save coords to file if required
              if(saveFileCoords.eq."Y")then
	            write(10,"(I8,I8,I8,F16.8,F16.8,F16.8)") i,m,innerKey(m),&
			    innerCoords(m,1),innerCoords(m,2),innerCoords(m,3)
			  endif
!store atom-force data
			  configurationRefForce(atomCounter,1) = configForcesR(n,1)
              configurationRefForce(atomCounter,2) = configForcesR(n,2)
              configurationRefForce(atomCounter,3) = configForcesR(n,3)			  
			  atomTypeKey(atomCounter) = configCoordsI(n)
			End Do
		  End Do
		End Do
      End Do
!----------------------------------------------------  
! x,y,z, min max
!----------------------------------------------------	  

!xMin
      vertexMinMax(1,1) = 0.0D0*workingUnitVector(1,1) + &
            0.0D0*workingUnitVector(1,2) + &
		    0.0D0*workingUnitVector(1,3)	!x1
      vertexMinMax(2,1) = 0.0D0*workingUnitVector(1,1) + &
            0.0D0*workingUnitVector(1,2) + &
	        1.0D0*alat*zCopy*workingUnitVector(1,3)	!x2
      vertexMinMax(3,1) = 0.0D0*workingUnitVector(1,1) + &
            1.0D0*alat*yCopy*workingUnitVector(1,2) + &
	        0.0D0*workingUnitVector(1,3)	!x3
      vertexMinMax(4,1) = 0.0D0*workingUnitVector(1,1) + &
            1.0D0*alat*yCopy*workingUnitVector(1,2) + &
	        1.0D0*alat*zCopy*workingUnitVector(1,3)	!x4
	  minMax(1) = vertexMinMax(1,1)
	  Do j=2,4		
		If(minMax(1).gt.vertexMinMax(j,1))Then
		  minMax(1) = vertexMinMax(j,1)
		End If
	  End Do	
	  minMax(1) = minMax(1)-1.5D0*rCutoff	  
!xMax
      vertexMinMax(5,1) = 1.0D0*alat*xCopy*workingUnitVector(1,1) + &
            0.0D0*workingUnitVector(1,2) + &
		    0.0D0*workingUnitVector(1,3)	!x1
      vertexMinMax(6,1) = 1.0D0*alat*xCopy*workingUnitVector(1,1) + &
            0.0D0*workingUnitVector(1,2) + &
	        1.0D0*alat*zCopy*workingUnitVector(1,3)	!x2
      vertexMinMax(7,1) = 1.0D0*alat*xCopy*workingUnitVector(1,1) + &
            1.0D0*alat*yCopy*workingUnitVector(1,2) + &
	        0.0D0*workingUnitVector(1,3)	!x3
      vertexMinMax(8,1) = 1.0D0*alat*xCopy*workingUnitVector(1,1) + &
            1.0D0*alat*yCopy*workingUnitVector(1,2) + &
	        1.0D0*alat*zCopy*workingUnitVector(1,3)	!x4
	  minMax(2) = vertexMinMax(5,1)
	  Do j=6,8		
		If(minMax(2).lt.vertexMinMax(j,1))Then
		  minMax(2) = vertexMinMax(j,1)
		End If
	  End Do	
	  minMax(2) = minMax(2)+1.5D0*rCutoff	
!yMin
      vertexMinMax(1,2) = 0.0D0*workingUnitVector(2,1) + &
            0.0D0*workingUnitVector(2,2) + &
		    0.0D0*workingUnitVector(2,3)	!x1
      vertexMinMax(2,2) = 0.0D0*workingUnitVector(2,1) + &
            0.0D0*workingUnitVector(2,2) + &
	        1.0D0*alat*zCopy*workingUnitVector(2,3)	!x2
      vertexMinMax(3,2) = 1.0D0*alat*xCopy*workingUnitVector(2,1) + &
            0.0D0*workingUnitVector(2,2) + &
	        0.0D0*workingUnitVector(2,3)	!x3
      vertexMinMax(4,2) = 1.0D0*alat*xCopy*workingUnitVector(2,1) + &
            0.0D0*workingUnitVector(2,2) + &
	        1.0D0*alat*zCopy*workingUnitVector(2,3)	!x4
	  minMax(3) = vertexMinMax(1,2)
	  Do j=2,4		
		If(minMax(3).gt.vertexMinMax(j,2))Then
		  minMax(3) = vertexMinMax(j,2)
		End If
	  End Do	
	  minMax(3) = minMax(3)-1.5D0*rCutoff	  
!yMax
      vertexMinMax(5,2) = 0.0D0*workingUnitVector(2,1) + &
            1.0D0*alat*yCopy*workingUnitVector(2,2) + &
		    0.0D0*workingUnitVector(2,3)	!x1
      vertexMinMax(6,2) = 0.0D0*workingUnitVector(2,1) + &
            1.0D0*alat*yCopy*workingUnitVector(2,2) + &
	        1.0D0*alat*zCopy*workingUnitVector(2,3)	!x2
      vertexMinMax(7,2) = 1.0D0*alat*xCopy*workingUnitVector(2,1) + &
            1.0D0*alat*yCopy*workingUnitVector(2,2) + &
	        0.0D0*workingUnitVector(2,3)	!x3
      vertexMinMax(8,2) = 1.0D0*alat*xCopy*workingUnitVector(2,1) + &
            1.0D0*alat*yCopy*workingUnitVector(2,2) + &
	        1.0D0*alat*zCopy*workingUnitVector(2,3)	!x4
	  minMax(4) = vertexMinMax(5,2)
	  Do j=6,8		
		If(minMax(4).lt.vertexMinMax(j,2))Then
		  minMax(4) = vertexMinMax(j,2)
		End If
	  End Do	
	  minMax(4) = minMax(4)+1.5D0*rCutoff		
!zMin
      vertexMinMax(1,3) = 0.0D0*workingUnitVector(3,1) + &
            0.0D0*workingUnitVector(3,2) + &
		    0.0D0*workingUnitVector(3,3)	!x1
      vertexMinMax(2,3) = 0.0D0*workingUnitVector(3,1) + &
            1.0D0*alat*yCopy*workingUnitVector(3,2) + &
	        0.0D0*workingUnitVector(3,3)	!x2
      vertexMinMax(3,3) = 1.0D0*alat*xCopy*workingUnitVector(3,1) + &
            0.0D0*workingUnitVector(3,2) + &
	        0.0D0*workingUnitVector(3,3)	!x3
      vertexMinMax(4,3) = 1.0D0*alat*xCopy*workingUnitVector(3,1) + &
            1.0D0*alat*yCopy*workingUnitVector(3,2) + &
	        0.0D0*workingUnitVector(3,3)	!x4   			
	  minMax(5) = vertexMinMax(1,3)
	  Do j=2,4		
		If(minMax(5).gt.vertexMinMax(j,3))Then
		  minMax(5) = vertexMinMax(j,3)
		End If
	  End Do	
	  minMax(5) = minMax(5)-1.5D0*rCutoff		
!zMax
      vertexMinMax(5,3) = 0.0D0*workingUnitVector(3,1) + &
            0.0D0*workingUnitVector(3,2) + &
		    1.0D0*alat*zCopy*workingUnitVector(3,3)	!x1
      vertexMinMax(6,3) = 1.0D0*alat*xCopy*workingUnitVector(3,1) + &
            0.0D0*workingUnitVector(3,2) + &
	        1.0D0*alat*zCopy*workingUnitVector(3,3)	!x2
      vertexMinMax(7,3) = 0.0D0*workingUnitVector(3,1) + &
            1.0D0*alat*yCopy*workingUnitVector(3,2) + &
	        1.0D0*alat*zCopy*workingUnitVector(3,3)	!x3
      vertexMinMax(8,3) = 1.0D0*alat*xCopy*workingUnitVector(3,1) + &
            1.0D0*alat*yCopy*workingUnitVector(3,2) + &
	        1.0D0*alat*zCopy*workingUnitVector(3,3)	!x4  
	  minMax(6) = vertexMinMax(5,3)
	  Do j=6,8		
		If(minMax(6).lt.vertexMinMax(j,3))Then
		  minMax(6) = vertexMinMax(j,3)
		End If
	  End Do	
	  minMax(6) = minMax(6)+1.5D0*rCutoff	

      print *, minMax(1), minMax(2)
      print *, minMax(3), minMax(4)
      print *, minMax(5), minMax(6)	  
	  
!----------------------------------------------------  
! Build Halo Volume
!---------------------------------------------------- 	
      j = 0
	  Do l=-1,1
	    Do m=-1,1
	      Do n=-1,1
	        If(l.eq.0.and.m.eq.0.and.n.eq.0)Then	
			  !Skip
			Else
	  mm = 0		
	  Do x=1,xCopy
	    Do y=1,yCopy
	      Do z=1,zCopy
		    Do k=coordsStart,(coordsStart+coordsLength-1)  
!calculate co-ordinates
			  xCoord = alat * (x + l*xCopy + configCoordsR(k,1) - 1)
			  yCoord = alat * (y + m*yCopy + configCoordsR(k,2) - 1)
			  zCoord = alat * (z + n*zCopy + configCoordsR(k,3) - 1)			  
			  !If(xCoord.ge.minMax(1).and.xCoord.le.minMax(2).and.&
			  !yCoord.ge.minMax(3).and.yCoord.le.minMax(4).and.&
			  !zCoord.ge.minMax(5).and.zCoord.le.minMax(6))Then	
			    mm = mm + 1	
			    j = j + 1	  
			    haloKey(j,1) = configCoordsI(k)		!atom type
			    haloKey(j,2) = mm					!real atom id
!x-coord
                haloCoords(j,1) = xCoord*workingUnitVector(1,1) + &
                yCoord*workingUnitVector(1,2) + zCoord*workingUnitVector(1,3)
!y-coord
                haloCoords(j,2) = xCoord * workingUnitVector(2,1) + &
                yCoord * workingUnitVector(2,2) + zCoord * workingUnitVector(2,3)
!z-coord
                haloCoords(j,3) = xCoord * workingUnitVector(3,1) + &
                yCoord * workingUnitVector(3,2) + zCoord * workingUnitVector(3,3)
			  !End If
			End Do
		  End Do
		End Do
	  End Do
			End If
		  End Do
	    End Do
	  End Do
	  haloAtoms = j
	  print *, atoms, haloAtoms
!----------------------------------------------------  
! Build Neighbour List
!---------------------------------------------------- 	 
      configLength = 0
!Compare with atoms in volume
      Do atomA=1,atoms
	    Do atomB=1,atoms
	      If(atomA.ne.atomB)Then
!check range one co-ord at a time then distance squared            
			If(atomA.gt.atomB)Then
			  atomAKey = atomB		!swap so atomA is lowest
			  atomBKey = atomA
			Else  
			  atomAKey = atomA
			  atomBKey = atomB
			End If
!atom coordinates
            xA = 1.0D0*innerCoords(atomAKey,1)
            xB = 1.0D0*innerCoords(atomBKey,1)
            xdSq = (xA-xB)**2
			If(xdSq.le.rCutoffSq)Then
              yA = 1.0D0*innerCoords(atomAKey,2)
              yB = 1.0D0*innerCoords(atomBKey,2)
			  ydSq = (yA-yB)**2
			  If(ydSq.le.rCutoffSq)Then
                zA = 1.0D0*innerCoords(atomAKey,3)
                zB = 1.0D0*innerCoords(atomBKey,3)
				zdSq = (zA-zB)**2
				If(zdSq.le.rCutoffSq)Then
				  rdSq = xdSq + ydSq + zdSq						
				  If(rdSq.le.rCutoffSq)Then
	                nlKey = nlKey + 1
					configLength = configLength + 1
					neighbourListITemp(nlKey,1) = innerKey(atomAKey)
					neighbourListITemp(nlKey,2) = innerKey(atomBKey)
					neighbourListITemp(nlKey,3) = atomAKey
					neighbourListITemp(nlKey,4) = atomBKey
					neighbourListCoordsTemp(nlKey,1) = xA
					neighbourListCoordsTemp(nlKey,2) = yA
					neighbourListCoordsTemp(nlKey,3) = zA
					neighbourListCoordsTemp(nlKey,4) = xB
					neighbourListCoordsTemp(nlKey,5) = yB
					neighbourListCoordsTemp(nlKey,6) = zB
					neighbourListRTemp(nlKey) = rdSq**0.5	
				  End If
				End If
	          End If
	        End If
	      End If
	    End Do
	  End Do	  
!Compare with atoms in halo
      Do atomA=1,atoms
	    Do atomB=1,haloAtoms
	      !If(atomA.ne.atomB)Then
!check range one co-ord at a time then distance squared   
			atomAKey = atomA		
			atomBKey = atoms+atomB
!atom coordinates
            xA = 1.0D0*innerCoords(atomAKey,1)
            xB = 1.0D0*haloCoords(atomB,1)
            xdSq = (xA-xB)**2
			If(xdSq.le.rCutoffSq)Then
              yA = 1.0D0*innerCoords(atomAKey,2)
              yB = 1.0D0*haloCoords(atomB,2)
			  ydSq = (yA-yB)**2
			  If(ydSq.le.rCutoffSq)Then
                zA = 1.0D0*innerCoords(atomAKey,3)
                zB = 1.0D0*haloCoords(atomB,3)
				zdSq = (zA-zB)**2
				If(zdSq.le.rCutoffSq)Then
				  rdSq = xdSq + ydSq + zdSq						
				  If(rdSq.le.rCutoffSq)Then
	                nlKey = nlKey + 1
					configLength = configLength + 1
					neighbourListITemp(nlKey,1) = innerKey(atomAKey)
					neighbourListITemp(nlKey,2) = haloKey(atomB,1)
					neighbourListITemp(nlKey,3) = atomAKey
					neighbourListITemp(nlKey,4) = haloKey(atomB,2)
					neighbourListCoordsTemp(nlKey,1) = xA
					neighbourListCoordsTemp(nlKey,2) = yA
					neighbourListCoordsTemp(nlKey,3) = zA
					neighbourListCoordsTemp(nlKey,4) = xB
					neighbourListCoordsTemp(nlKey,5) = yB
					neighbourListCoordsTemp(nlKey,6) = zB
					neighbourListRTemp(nlKey) = rdSq**0.5	
				  End If
				End If
	          End If
	        End If
	      !End If
	    End Do
	  End Do	