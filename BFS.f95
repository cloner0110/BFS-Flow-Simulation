PROGRAM BACKWARD
IMPLICIT NONE
INTEGER :: I,J,K,GRID1,GRID2,GRID_TOTAL,ITER,O,ITERATION
INTEGER :: NX1,NY1,NY2,NX2,NY,NX
REAL :: DX,DY,L1,L2,H1,H2,ULEFT,VLEFT,L,H,XP,YP
REAL :: EPSILON_C,NU,RO,EPSILON_P
REAL :: P_P1,P_P2,P_P3,P_P4,P_P5,P_P6
REAL :: P_NX1,P_NX2,P_NX3,P_NX4,P_NX5,P_NX6,P_NX7
REAL :: P_NY1,P_NY2,P_NY3,P_NY4,P_NY5,P_NY6,P_NY7
REAL , ALLOCATABLE :: U(:,:)
REAL , ALLOCATABLE :: V(:,:)
REAL , ALLOCATABLE :: P(:,:)
REAL , ALLOCATABLE :: UTEMP(:,:)
REAL , ALLOCATABLE :: VTEMP(:,:)
REAL , ALLOCATABLE :: UINIT(:,:)
REAL , ALLOCATABLE :: VINIT(:,:)
REAL , ALLOCATABLE :: PTEMP(:,:)
REAL , ALLOCATABLE :: POLD(:,:)
REAL , ALLOCATABLE :: PNEW(:,:)
REAL , ALLOCATABLE :: RE(:,:)
REAL , ALLOCATABLE :: UV(:,:)
REAL , ALLOCATABLE :: UP(:,:)
REAL , ALLOCATABLE :: VP(:,:)
REAL , ALLOCATABLE :: D(:,:)
REAL , ALLOCATABLE :: X(:,:)
REAL , ALLOCATABLE :: Y(:,:)
!**************************
!H=H1+H2
!L=L1+L2
!NX=NX1+NX2-1
!NY=NY1+NY2-1
!*************************************
!GRID1=NY1*(NX1+(NX2-1))
!GRID2=NX2*(NY2-1)
!GRID_TOTAL=GRID1+GRID2
!=========================================
RO=998.2
NU=0.001003/RO
DX=0.9
DY=0.1
NX1=20
NX2=80
NY1=7
NY2=10
ULEFT=0.00001
L1=NX1*DX
L2=(NX2-1)*DX
H1=NY1*DY
H2=(NY2-1)*DY
EPSILON_C=0.01
EPSILON_P=0.01
H=H1+H2
L=L1+L2
NX=NX1+NX2-1
NY=NY1+NY2-1
!=========================================
allocate (U(NX,NY))
allocate (UTEMP(NX,NY))
allocate (P(NX,NY))
allocate (V(NX,NY))
allocate (VTEMP(NX,NY))
allocate (PTEMP(NX,NY))
allocate (UV(NX,NY))
allocate (UP(NX,NY))
allocate (VP(NX,NY))
allocate (D(NX,NY))
allocate (UINIT(NX,NY))
allocate (VINIT(NX,NY))
allocate (X(NX,NY))
allocate (Y(NX,NY))
allocate (RE(NX,NY))
allocate (POLD(NX,NY))
allocate (PNEW(NX,NY))
!************************** INITIAL BOUNDARY VALUE 
UINIT(:,:)=0.5
VINIT(:,:)=0.5
DO I=1,NX
	DO J=1,NY
	X(I,J)=(I-1)*DX
	Y(I,J)=(J-1)*DY
	END DO
END DO
DO I=1,NX
	DO J=1,NY
		IF ( X(I,J)>=0 .AND. X(I,J)<=L1 .AND. Y(I,J)>=0 .AND. Y(I,J)<=H1 ) THEN
		U(I,J)=0.
		V(I,J)=0.
		END IF
		IF (X(I,J)==0 .AND. Y(I,J)>H1 .AND. Y(I,J)<=H) THEN
		U(1,J)=ULEFT
		V(1,J)=0
		END IF 
		IF (Y(I,J)==0 .AND. X(I,J)>=L1 .AND. X(I,J)<=L ) THEN
		U(I,J)=0
		V(I,J)=0
		END IF
	END DO
END DO
DO I=1,NX
  U(I,NY)=0
  V(I,NY)=0
END DO
!===================================== .REYNOLDS
DO I=1,NX
  DO J=1,NY
    IF ( X(I,J)<=L1 ) THEN
      RE(I,J)=(ULEFT*H2)/NU
    ELSE 
      RE(I,J)=((H2/H)**2)*ULEFT
    END IF
  END DO
END DO
!===================================== ..EQUATIONS..
!UINIT(:,:)=0
!VINIT(:,:)=0
PTEMP(:,:)=0
UTEMP(:,:)=UINIT(:,:)
VTEMP(:,:)=VINIT(:,:)
 DO ITER=1,150000
  DO I=2,NX-1
    DO J=2,NY-1
		IF(J<=NY1 .AND. I<=NX1 ) THEN
        		GOTO 10
        	END IF
P_NY1=(PTEMP(I,J)-PTEMP(I,J-1))/(RO*DY)
P_NY2=(UTEMP(I,J)*(VTEMP(I,J)-VTEMP(I-1,J)))/DX
P_NY3=(VTEMP(I,J)-VTEMP(I,J-1))/DY
P_NY4=VTEMP(I,J)/DY
P_NY5=(VTEMP(I,J)*VTEMP(I,J-1))/DY
P_NY6=2*DY*(VTEMP(I+1,J)+VTEMP(I-1,J))
P_NY7=2*DX*(VTEMP(I,J+1)+VTEMP(I,J-1))
!******************************************************
P_NX1=(PTEMP(I,J)-PTEMP(I-1,J))/(RO*DX)
P_NX2=(VTEMP(I,J)*(UTEMP(I,J)-UTEMP(I,J-1)))/DY
P_NX3=(UTEMP(I,J)-UTEMP(I-1,J))/DX
P_NX4=UTEMP(I,J)/DX
P_NX5=(UTEMP(I,J)*UTEMP(I-1,J))/DX
P_NX6=2*DY*(UTEMP(I+1,J)+UTEMP(I-1,J))
P_NX7=2*DX*(UTEMP(I,J+1)+UTEMP(I,J-1))
!******************************************************		
U(I,J)=((P_NX5)-(P_NX1)-(P_NX2)+(NU/4*DX*DY)*(P_NX6+P_NX7))/((P_NX3)+(P_NX4)+(NU*(DX+DY))/(DX*DY))
V(I,J)=((P_NY5)-(P_NY1)-(P_NY2)+(NU/4*DX*DY)*(P_NY6+P_NY7))/((P_NY3)+(P_NY4)+(NU*(DX+DY))/(DX*DY))
10  END DO
END DO
!============================= SHARTE MARZI DAR KHORUJI , GRADIAN SORAT =0
DO J=1,NY
  U(NX,J)=U(NX-1,J)
  V(NX,J)=V(NX-1,J)
  END DO
!=======================================   
D(:,:)=0
UP(:,:)=0
VP(:,:)=0
UV(:,:)=0       
   DO I=2,NX-1
     DO J=2,NY-1
       D(I,J)=((U(I,J)-U(I-1,J))/DX)+((V(I,J)-V(I,J-1))/DY) 
      END DO
      END DO
DO I=2,NX-1
  DO J=2,NY-1
    UP(I,J)=U(I,J)**2
    VP(I,J)=V(I,J)**2
    UV(I,J)=U(I,J)*V(I,J)
    END DO
    END DO
 !======================================== solving pressure
POLD(:,:)=0
PNEW(:,:)=0    
DO ITERATION=1,1500000
DO I=2,NX-1
  DO J=2,NY-1
P_P1=2*(DX+DY)
P_P2=DY*(POLD(I+1,J)+POLD(I-1,J))+DX*(POLD(I,J+1)+POLD(I,J-1))
P_P3=(D(I+1,J)-2*D(I,J)+D(I-1,J))/(2*DX) +(D(I,J+1)-2*D(I,J)+D(I,J-1))/(2*DY)
P_P4=(UP(I+1,J)-2*UP(I,J)+UP(I-1,J))/(2*DX)
P_P5=0.5*(UV(I+1,J+1)-UV(I-1,J+1)-UV(I+1,J-1)+UV(I-1,J-1)) 
P_P6=(VP(I+1,J)-2*VP(I,J)+VP(I-1,J))/(2*DY) 
PNEW(I,J)= (P_P2-(2*DX*DY)*(P_P3/RE(I,J))-(P_P4)-(P_P5)-(P_P6))/P_P1
! PRINT*,PNEW(I,J),"FESHAR"
END DO
END DO
O=0
DO I=2,NX-1
  DO J=2,NY-1
   IF (ABS(PNEW(I,J)-POLD(I,J))<=EPSILON_P ) THEN
     O=O+1  
END IF
END DO
END DO
POLD(:,:)=PNEW(:,:)                      
if (o==NX*NY)THEN
   GOTO 100                      ! DOROSTE ?
  END IF
  END DO
100 P(:,:)=PNEW(:,:)
!============================================= ! checking the Continuity with D
O=0
  DO I=2,NX-1    ! SOAL
  DO J=2,NY-1
  IF (((U(I+1,J)-2*U(I,J)+U(I-1,J))/(2*DX)+(V(I,J+1)-2*V(I,J)+V(I,J-1))/(2*DY))<EPSILON_C) THEN
	O=O+1
    PRINT*,o
END IF		
		END DO 
        END DO
IF ( O==(NX)*(NY) ) THEN
  GOTO 13        
     END IF
     	UTEMP(:,:)=U(:,:)
      	VTEMP(:,:)=V(:,:)
      	PTEMP(:,:)=P(:,:) 
     END DO         					
   13   UTEMP(:,:)=U(:,:)
      	VTEMP(:,:)=V(:,:)
      	PTEMP(:,:)=P(:,:) 

     
   DO ITER=1,150000
 	 DO I=2,NX-1
  14  DO J=2,NY-1
		IF(J<=NY1 .AND. I<=NX1 ) THEN
        GOTO 14
        END IF
P_NY1=(PTEMP(I,J)-PTEMP(I,J-1))/(RO*DY)
P_NY2=(UTEMP(I,J)*(VTEMP(I,J)-VTEMP(I-1,J)))/DX
P_NY3=(VTEMP(I,J)-VTEMP(I,J-1))/DY
P_NY4=VTEMP(I,J)/DY
P_NY5=(VTEMP(I,J)*VTEMP(I,J-1))/DY
P_NY6=2*DY*(VTEMP(I+1,J)+VTEMP(I-1,J))
P_NY7=2*DX*(VTEMP(I,J+1)+VTEMP(I,J-1))
!******************************************************
P_NX1=(PTEMP(I,J)-PTEMP(I-1,J))/(RO*DX)
P_NX2=(VTEMP(I,J)*(UTEMP(I,J)-UTEMP(I,J-1)))/DY
P_NX3=(UTEMP(I,J)-UTEMP(I-1,J))/DX
P_NX4=UTEMP(I,J)/DX
P_NX5=(UTEMP(I,J)*UTEMP(I-1,J))/DX
P_NY6=2*DY*(UTEMP(I+1,J)+UTEMP(I-1,J))
P_NY7=2*DX*(UTEMP(I,J+1)+UTEMP(I,J-1))
!******************************************************		
U(I,J)=((P_NX5)-(P_NX1)-(P_NX2)+(NU/4*DX*DY)*(P_NX6+P_NX7))/((P_NX3)+(P_NX4)+(NU*(DX+DY))/(DX*DY))
V(I,J)=((P_NY5)-(P_NY1)-(P_NY2)+(NU/4*DX*DY)*(P_NY6+P_NY7))/((P_NY3)+(P_NY4)+(NU*(DX+DY))/(DX*DY))
        END DO
        END DO    
DO I=2,NX-1
  DO J=2,NY-1
    UP(I,J)=U(I,J)**2
    VP(I,J)=V(I,J)**2
    UV(I,J)=U(I,J)*V(I,J)
    END DO
    END DO
O=0
POLD(:,:)=0
DO ITERATION=1,1500000  
DO I=1,NX
  DO J=1,NY
P_P1=2*(DX+DY)
P_P2=DY*(POLD(I+1,J)+POLD(I-1,J))+DX*(POLD(I,J+1)+POLD(I,J-1))
P_P3=(D(I+1,J)-2*D(I,J)+D(I-1,J))/(2*DX) +(D(I,J+1)-2*D(I,J)+D(I,J-1))/(2*DY)
P_P4=(UP(I+1,J)-2*UP(I,J)+UP(I-1,J))/(2*DX)
P_P5=0.5*(UV(I+1,J+1)-UV(I-1,J+1)-UV(I+1,J-1)+UV(I-1,J-1)) 
P_P6=(VP(I+1,J)-2*VP(I,J)+VP(I-1,J))/(2*DY) 
PNEW(I,J)= (P_P2-(2*DX*DY)*(P_P3/RE(I,J))-(P_P4)-(P_P5)-(P_P6))/P_P1  
END DO
END DO
DO I=1,NX
  DO J=1,NY
   IF (ABS(PNEW(I,J)-POLD(I,J))<=EPSILON_P ) THEN
     O=O+1
     PRINT*,"SALAM"
END IF
END DO
END DO
POLD(:,:)=PNEW(:,:)                     
if (o==NX*NY) exit                     
  END DO
P(:,:)=PNEW(:,:)
OPEN(40,FILE='VELOCITY.TXT')
WRITE(40,*)'VARIABLES=,"X","Y","U","V"'
WRITE(40,*)'ZONE I=','NX','J=','NY','F=POINT'      
 DO J=1,NY
        DO I=1,NX
      WRITE(40,200)X(I,J),Y(I,J),U(I,J),V(I,J)
      200 FORMAT(4F25.16)
      END DO
      END DO
END PROGRAM BACKWARD   



 
