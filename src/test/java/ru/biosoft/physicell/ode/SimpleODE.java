package ru.biosoft.physicell.ode;


import java.util.Map;

public class SimpleODE extends ODEModel
{
    public double rate_Reaction_1;
    public double RATE_OF_X;
    public double assignment;
    public double k;
    public double[] getY()
    {
        return x_values;
    }


    private void calculateParameters() throws Exception
    {
        double[] x_values = this.x_values;
        RATE_OF_X = -rate_Reaction_1;
    }


    private void calculateReactionRates()
    {
        double[] x_values = this.x_values;
        rate_Reaction_1 = k;
    }


    private void calculateInitialParameters()
    {
        double[] x_values = this.x_values;
        rate_Reaction_1 = k;
        RATE_OF_X = -rate_Reaction_1;
    }


    public final double[] dy_dt_slow(double time, double[] x_values) throws Exception
    {
        this.time = time;
        this.x_values = x_values;
        final double[] dydt = new double[1];
        calculateParameters();
        calculateReactionRates();
        dydt[0] = -rate_Reaction_1; //  rate rule for of $X
        return dydt;
    }




    @Override
    public final void init() throws Exception
    {
        CONSTRAINTS__VIOLATED = 0;
        this.simulationResultHistory.clear();
        this.simulationResultTimes.clear();
        rate_Reaction_1 = 0.0; // initial value of $$rate_Reaction_1
        RATE_OF_X = 0.0; // initial value of RATE_OF_X
        assignment = 0.0; // initial value of assignment
        k = 1.0; // initial value of k
        time = 0.0; // initial value of time
        calculateInitialValues();
        initMap();
        this.isInit = true;
    }


    @Override
    public final void init(double[] initialValues, Map parameters) throws Exception
    {
        super.init( initialValues, parameters );
        this.initialValues = x_values.clone();
    }


    private final void calculateInitialValues() throws Exception
    {
        double[] x_values = this.x_values = new double[1];
        this.time = 0.0;
        x_values[0] = 10.0; //  initial value of $X
        this.initialValues = x_values;
        calculateInitialParameters();
        this.initialValues = x_values;
    }




    public final double[] extendResult(double time, double[] x_values) throws Exception
    {
        this.time = time;
        this.x_values = x_values;
        calculateParameters();
        return getCurrentState();
    }


    public final double[] getCurrentState()
    {
        return new double[] {rate_Reaction_1, x_values[0], k, time,};
    }


    @Override
    public final void setCurrentValues(double[] values) throws Exception
    {
        CONSTRAINTS__VIOLATED = 0;
        rate_Reaction_1 = values[0];
        x_values[0] = values[1];
        k = values[2];
        time = values[3];
        if( time == 0 )
        {
            initialValues[0] = values[1];
            calculateInitialParameters();
        }
        else
            calculateParameters();
    }




    public final double[] checkEvent(double time, double[] x_values) throws Exception
    {
        this.time = time;
        this.x_values = x_values;
        calculateParameters();
        double[] z = new double[2];
        z[0] = ( NUMERIC_LEQ( x_values[0], 1 ) != 0 ) ? 1 : -1; // event_1 part1;
        z[1] = ( NUMERIC_LEQ( x_values[0], 1 ) != 0 ) ? 1 : -1; // event_1 part2;
        return z;
    }


    public final void processEvent(int i)
    {
        double[] assignments;
        double executionTime;
        switch( i )
        {
            case ( 0 ): //event_1 part1
            {
                assignments = new double[2 - 1];
                executionTime = time;
                assignments[1 - 1] = 10.0;
                addDelayedEvent( 0 + 1, executionTime, assignments );
                break;
            }
            case ( 1 ): //event_1 part2
            {
                assignments = getNextAssignments( 1 );
                x_values[0] = assignments[0] * 1.0;
                removeDelayedEvent( 1 );
                break;
            }
            default:
                break;
        }
    }


    public final double[] getEventsPriority(double time, double[] x_values) throws Exception
    {
        calculateParameters();
        return new double[] {Double.POSITIVE_INFINITY, //event_1 part1
                Double.NEGATIVE_INFINITY, //event_1 part2
        };
    }


    public final boolean getEventsInitialValue(int i) throws IndexOutOfBoundsException
    {
        return true;
    }


    public final boolean isEventTriggerPersistent(int i) throws IndexOutOfBoundsException
    {
        return true;
    }


    public final String getEventMessage(int i) throws IndexOutOfBoundsException
    {
        return null;
    }

    @Override
    protected void initMap()
    {
        variableIndex.put( "X", 1 );
    }

}
