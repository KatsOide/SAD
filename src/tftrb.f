      subroutine tftrb(latt,l1,sv)
      use ffs
      use tffitcode
      dimension sv(5),latt(2,l1)
      call tracke(latt,l1,sv,np,'STANDBY',' ',0)
      call tracke(latt,l1,sv,np,'TRACK',' ',0)
      call tracke(latt,l1,sv,np,'RESET',' ',0)
      return
      end
