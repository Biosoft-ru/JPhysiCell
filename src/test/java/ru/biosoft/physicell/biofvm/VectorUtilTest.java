package ru.biosoft.physicell.biofvm;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import org.junit.Test;

public class VectorUtilTest
{

    @Test
    public void testDouble()
    {
        double[] arr = null; //Create null array
        arr = VectorUtil.resize( arr, 5 ); //Resize array to length 5
        assertNotNull( arr );
        assertEquals( arr.length, 5 );

        for( int i = 0; i < arr.length; i++ ) //Initialize array
            arr[i] = i;

        arr = VectorUtil.resize( arr, 4 ); //Decrease array length
        assertEquals( arr.length, 4 );
        assert ( arr[3] == 3.0 );

        arr = VectorUtil.resize( arr, 6, 89 ); //Increase array length with given value
        assert ( arr[5] == 89 );

        arr = VectorUtil.resize( arr, 8 ); //Increase array length
        assert ( arr[7] == 0 );

        arr = VectorUtil.assign( 10, 1.33 ); //Initialize via assign
        assertEquals( arr.length, 10 );
        for( int i = 0; i < arr.length; i++ )
            assert ( arr[i] == 1.33 );

        //push back test
        arr = null;
        arr = VectorUtil.push_back( arr, 15 );
        assertEquals( arr.length, 1 );
        assert ( arr[0] == 15 );
        arr = VectorUtil.push_back( arr, 35 );
        assertEquals( arr.length, 2 );
        assert ( arr[1] == 35 );

        //equals test
        assert ( !VectorUtil.equals( new double[] {2.2, 2.1}, new double[] {2.1, 2.2} ) );
        assert ( !VectorUtil.equals( new double[] {2.2, 2.1}, new double[] {2.2} ) );
        assert ( !VectorUtil.equals( new double[] {2.2}, new double[] {2.2, 10} ) );
        assertEquals( new double[] {2.2, 10}, new double[] {2.2, 10} );

        //array of arrays
        arr = new double[] {3, 4};
        double[][] arrs = VectorUtil.assign( 3, arr );
        arrs[0][0] = 6;
        assertEquals( arrs.length, 3 );
        assert ( arrs[1][0] == 3 ); //check that it was not affected by arr[0][0] setting

        arrs = VectorUtil.resize( arrs, 4 );// increase array length
        assertEquals( arrs.length, 4 );
        assertEquals( arrs[3].length, 2 );
        assert ( arrs[1][0] == 3 );

        arrs = VectorUtil.resize( arrs, 2 ); //decrease array length
        assertEquals( arrs.length, 2 );
        assert ( arrs[1][0] == 3 );
    }

    @Test
    public void testBoolean()
    {
        boolean[] arr = null; //Create null array
        arr = VectorUtil.resize( arr, 5 ); //Resize array to length 5
        assertNotNull( arr );
        assertEquals( arr.length, 5 );

        for( int i = 0; i < arr.length; i++ ) //Initialize array
            arr[i] = true;

        arr = VectorUtil.resize( arr, 4 ); //Decrease array length
        assertEquals( arr.length, 4 );
        assert ( arr[3] );

        arr = VectorUtil.resize( arr, 6 ); //Increase array length
        assert ( !arr[5] );

        //push back test
        arr = null;
        arr = VectorUtil.push_back( arr, true );
        assertEquals( arr.length, 1 );
        assert ( arr[0] );
        arr = VectorUtil.push_back( arr, false );
        assertEquals( arr.length, 2 );
        assert ( !arr[1] );

        //array of arrays
        arr = new boolean[] {true, true};
        boolean[][] arrs = VectorUtil.assign( 3, arr );
        arrs[0][0] = false;
        assertEquals( arrs.length, 3 );
        assert ( arrs[1][0] ); //check that it was not affected by arr[0][0] setting
    }

    @Test
    public void testArithemtic()
    {
        //sum
        double[] arr1 = VectorUtil.assign( 2, 1.1 );
        double[] arr2 = new double[2];

        VectorUtil.sum( arr2, arr1 ); //arr2 = arr2 + arr1
        assertEquals( arr2, arr1 );
        VectorUtil.sum( arr2, new double[] {1.1, 1.1} ); //arr2 = arr2 + arr1 + arr1
        assertEquals( arr2, new double[] {2.2, 2.2} );

        assertEquals( VectorUtil.newSum( arr1, 1.2 ), new double[] {2.3, 2.3} ); //arr1 + 1.2
        assertEquals( VectorUtil.newSum( 1.2, arr1 ), new double[] {2.3, 2.3} ); //1.3 + arr1
        assertEquals( VectorUtil.newSum( arr1, new double[] {1.0, 4.0} ), new double[] {2.1, 5.1} ); //arr1 + array
        assertEquals( VectorUtil.newSum( arr1, arr2 ), VectorUtil.newSum( arr2, arr1 ) ); //symmetric

        //product
        arr1 = new double[] {1, 1};
        VectorUtil.prod( arr1, 2 ); //arr1 *= 2
        assertEuals( arr1, new double[] {2, 2} );
        assertEuals( VectorUtil.newProd( arr1, 2 ), new double[] {4, 4} ); //arr1 * 2
        assertEuals( VectorUtil.newProd( 3, arr1 ), new double[] {6, 6} ); //3 * arr1
        assertEuals( VectorUtil.newProd( arr1, new double[] {2, 4} ), new double[] {4, 8} ); // arr1 * [2,4] = [4,8]

        //difference
        arr1 = new double[] {3, 3};
        arr2 = new double[] {1, 0};
        VectorUtil.diff( 5, arr1 );
        assertEuals( arr1, new double[] {2, 2} );
        VectorUtil.diff( arr1, arr2 );
        assertEuals( arr1, new double[] {1, 2} );

        assertEuals( VectorUtil.newDiff( arr1, 1 ), new double[] {0, 1} );
        assertEuals( VectorUtil.newDiff( arr1, new double[] {0.5, 0.1} ), new double[] {0.5, 1.9} );

        //division
        arr1 = new double[] {2, 4};
        arr2 = new double[] {16, 8};
        VectorUtil.div( arr2, 2 );
        assertEuals( arr2, new double[] {8, 4} );
        VectorUtil.div( arr2, arr1 );
        assertEuals( arr2, new double[] {4, 1} );
        assertEuals( VectorUtil.newDiv( arr2, arr1 ), new double[] {2, 0.25} );

        //BLAS
        arr1 = new double[] {5, 6};
        arr2 = new double[] {0, 2};
        double[] arr3 = new double[] {2, 3};
        VectorUtil.axpy( arr1, 2, arr2 );
        assertEuals( arr1, new double[] {5, 10} );
        VectorUtil.naxpy( arr1, 3, arr2 );
        assertEuals( arr1, new double[] {5, 4} );
        VectorUtil.axpy( arr1, arr2, arr3 );
        assertEuals( arr1, new double[] {5, 10} );
        VectorUtil.naxpy( arr1, arr3, arr3 );
        assertEuals( arr1, new double[] {1, 1} );
    }

    private void assertEuals(double[] arr1, double[] arr2)
    {
        assert ( VectorUtil.equals( arr1, arr2 ) );
    }
}
