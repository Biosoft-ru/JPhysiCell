package ru.biosoft.physicell.core;

import org.junit.Test;

import ru.biosoft.physicell.biofvm.Microenvironment;

public class MicroenvironmentTest
{


    @Test
    public void testSimple()
    {
        //Create simple microenvironment with uniform mesh (5x5x5)
        Microenvironment m = new Microenvironment( "substrate scale", 100, 20, "minutes", "microns" );
        m.setDensity( 0, "oxygen", "mmHg", 0, 0 );
        assert ( m.numberVoxels() == 125 );

        assert ( m.get( 0 ).length == 1 ); //at ith coordinate vector of length 1 (oxygen)

        //Add more density and resize to more complex mesh (20x10x5)
        m.addDensity( "oxygen2", "mmHGg", 10, 10 );
        m.resizeSpace( -100, 100, 0, 10, 100, 150, 10, 1, 10 );
        assert ( m.numberVoxels() == 1000 );
        assert ( m.get( 0 ).length == 2 );
    }
}
