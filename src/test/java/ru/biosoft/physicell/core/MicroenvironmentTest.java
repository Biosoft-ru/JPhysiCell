package ru.biosoft.physicell.core;

import org.junit.Test;

import ru.biosoft.physicell.biofvm.Microenvironment;

public class MicroenvironmentTest
{

    /**
     * Create simple environment with 5x5x5 mesh 
     */
    @Test
    public void test()
    {
        Microenvironment m = new Microenvironment( "substrate scale", 100, 20, "minutes", "microns" );
        m.setDensity( 0, "oxygen", "mmHg" );
        assert ( m.number_of_voxels() == 125 );
    }
}