package ru.biosoft.physicell.ode;

import biouml.plugins.simulation.java.JavaBaseModel;

public class Euler
{
    protected double t; //current model time value
    protected double[] xOld; //vlaues of ODE variables at previous step
    protected double[] x; //current values of all ODE variables
    protected double tFinal; //final time
    protected double h; // time step (variable)
    protected double hBase;
    protected double nextTime;
    protected double dt;
    boolean eventOccured = false;

    JavaBaseModel model;

    public Euler(ODEModel model) throws Exception
    {
        this( model, 0.01 );
    }

    public Euler(ODEModel model, double dt) throws Exception
    {
        this.model = model;
        if( !model.isInit() )
            model.init();

        x = model.getInitialValues();
        this.dt = dt;
        this.h = dt;
        this.hBase = dt;
        this.nextTime = dt;
        prevEvent = model.checkEvent( t, x );
        t = 0;//span.getTimeStart();
    }

    double[] prevEvent;

    public void doStep() throws Exception
    {
        if (t == 0)
            initialEvent();
        while( t < nextTime )
        {
            xOld = model.x_values.clone();
            h = Math.min( hBase, nextTime - t ); //restrict step size, so we won't get beyond span point
            integrationStep( x, xOld, t, h );
            t += h;
            model.extendResult( t, x );

            double[] ev = model.checkEvent( t, x );
            eventOccured = false;
            for( int i = 0; i < ev.length; i++ )
            {
                if( ev[i] > 0 && prevEvent[i] < 0 )
                {
                    eventOccured = true;
                    model.processEvent( i );
                }
            }
            if( eventOccured )
                model.extendResult( t, x );
            ev = prevEvent.clone();
        }
        nextTime += dt;
    }

    private void initialEvent() throws Exception
    {
        double[] prevEvent = model.checkEvent( 0, x );
        for( int i = 0; i < prevEvent.length; i++ )
        {
            if( prevEvent[i] > 0 )
            {
                eventOccured = true;
                //                System.out.println( "Event " + i );
                model.processEvent( i );
            }
        }
        if( eventOccured )
            model.extendResult( t, x );
    }

    public void integrationStep(double[] xNew, double[] xOld, double tOld, double h, double theta) throws Exception
    {
        integrationStep( xNew, xOld, tOld, h * theta );
    }

    public void integrationStep(double[] xNew, double[] xOld, double tOld, double h) throws Exception
    {
        double[] dydt = model.dy_dt( tOld, xOld );

        for( int i = 0; i < xNew.length; i++ )
            xNew[i] = xOld[i] + dydt[i] * h;
    }

    public Euler clone(ODEModel model) throws Exception
    {
        Euler result = new Euler( model );
        result.dt = dt;
        return result;
    }
}