import itk
import sys

if len (sys.argv) < 3:
    print( "Usage: %s <input> <output>" % (sys.argv[0]) )
    sys.exit(1)

image = itk.imread(sys.argv[1], itk.UC)
filled = itk.morphological_contour_interpolator(image)
itk.imwrite(filled, sys.argv[2], compression=True)
