package ru.biosoft.physicell.ode;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class TestSimpleODE
{
    @Test
    public void test() throws Exception
    {
        SimpleODE model = new SimpleODE();
        Euler euler = new Euler( model, 0.1 );
        double initial = 10;
        assert ( model.getCurrentValue( "X" ) == initial );

        for( int i = 0; i < 90; i++ )
        {
            euler.doStep();
            System.out.println( model.getCurrentValue( "X" ) );

        }

        euler.doStep();
        assertEquals( model.getCurrentValue( "X" ), initial, 1E-5 );
    }
}
