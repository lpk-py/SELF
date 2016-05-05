PROGRAM BSW_Driver


 USE ModelPrecision
 USE ModelFlags
 USE ConstantsDictionary

 USE BarotropicShelfWaves_Class

 IMPLICIT NONE

 TYPE( BarotropicShelfWaves ) :: mybsw

    CALL mybsw % Build( )

    CALL mybsw % IRAM( )

    CALL mybsw % WriteTecplot( )

    CALL mybsw % Trash( )

END PROGRAM BSW_Driver
