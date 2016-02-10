      subroutine tins(np,x,px,y,py,z,g,trans)
      include 'inc/TMACRO.inc'
      real*8 x(np),px(np),y(np),py(np),z(np),g(np)
      real*8 trans(6,7)
      do i=1,np
c        dp=g(i)*(2.d0+g(i))
        dp=g(i)
        pr=1.d0+dp
        x1=x(i)
        px1=px(i)*pr
        y1=y(i)
        py1=py(i)*pr
        x(i) = trans(1,1)*x1+trans(1,2)*px1+
     $         trans(1,3)*y1+trans(1,4)*py1+
     $         trans(1,6)*dp+trans(1,7)
        px(i)=(trans(2,1)*x1+trans(2,2)*px1+
     $         trans(2,3)*y1+trans(2,4)*py1+
     $         trans(2,6)*dp+trans(2,7))/pr
        y(i) = trans(3,1)*x1+trans(3,2)*px1+
     $         trans(3,3)*y1+trans(3,4)*py1+
     $         trans(3,6)*dp+trans(3,7)
        py(i)=(trans(4,1)*x1+trans(4,2)*px1+
     $         trans(4,3)*y1+trans(4,4)*py1+
     $         trans(4,6)*dp+trans(4,7))/pr
        z(i) = trans(5,1)*x1+trans(5,2)*px1+
     $         trans(5,3)*y1+trans(5,4)*py1+trans(5,6)*dp+z(i)
      enddo
      return
      end
