package ru.biosoft.physicell.sample_projects.ode_energy;

/*
* This code is generated by BioUML FrameWork
* for Toy_Metabolic_Model diagram  at $creationTime
*/
import java.util.Map;

import ru.biosoft.math.MathRoutines;
import ru.biosoft.physicell.ode.ODEModel;

public class ToyMetabolicModel extends ODEModel
{
    public double rate_Aerobic;
    public double rate_Anaerobic;
    public double rate_Energy_Usage;
    public double Intracellular;
    public double Intracellular_Lac_Secretion_Rate;
    public double Intracellular_Transition_Rate;
    public double Intracellular_apoptosis_rate;
    public double Intracellular_migration_speed;
    public double RATE_OF_Intracellular_Energy;
    public double RATE_OF_Intracellular_Glucose;
    public double RATE_OF_Intracellular_Lactate;
    public double RATE_OF_Intracellular_Oxygen;
    public double assignment;
    public double energy_death_thresh;
    public double energy_prolif_thresh;
    public double k_aer;
    public double k_ane;
    public double k_usage;
    public double oxygen_thresh;
    public double[] getY()
    {
        return x_values;
    }


    private void calculateParameters() throws Exception
    {
        double[] x_values = this.x_values;
        RATE_OF_Intracellular_Glucose = ( -rate_Aerobic - rate_Anaerobic ) / Intracellular;
        RATE_OF_Intracellular_Energy = ( rate_Aerobic * 38 + rate_Anaerobic * 2 - rate_Energy_Usage ) / Intracellular;
        RATE_OF_Intracellular_Lactate = rate_Anaerobic / Intracellular;
        RATE_OF_Intracellular_Oxygen = -rate_Aerobic * 6 / Intracellular;
    }


    private void calculateReactionRates()
    {
        double[] x_values = this.x_values;
        rate_Aerobic = Intracellular * k_aer * ( x_values[1] / Intracellular )
                * MathRoutines.pow( ( x_values[3] / Intracellular ), (int)6.0 );
        rate_Anaerobic = Intracellular * k_ane * ( x_values[1] / Intracellular );
        rate_Energy_Usage = Intracellular * k_usage * ( x_values[0] / Intracellular );
    }


    private void calculateInitialParameters()
    {
        double[] x_values = this.x_values;
        x_values[1] = initialValues[1] * Intracellular;
        x_values[3] = initialValues[3] * Intracellular;
        Intracellular_Lac_Secretion_Rate = Intracellular_Lac_Secretion_Rate * Intracellular;
        Intracellular_Transition_Rate = Intracellular_Transition_Rate * Intracellular;
        x_values[2] = initialValues[2] * Intracellular;
        x_values[0] = initialValues[0] * Intracellular;
        Intracellular_migration_speed = Intracellular_migration_speed * Intracellular;
        Intracellular_apoptosis_rate = Intracellular_apoptosis_rate * Intracellular;
        rate_Aerobic = Intracellular * k_aer * ( x_values[1] / Intracellular )
                * MathRoutines.pow( ( x_values[3] / Intracellular ), (int)6.0 );
        RATE_OF_Intracellular_Oxygen = -rate_Aerobic * 6 / Intracellular;
        rate_Anaerobic = Intracellular * k_ane * ( x_values[1] / Intracellular );
        RATE_OF_Intracellular_Lactate = rate_Anaerobic / Intracellular;
        RATE_OF_Intracellular_Glucose = ( -rate_Aerobic - rate_Anaerobic ) / Intracellular;
        rate_Energy_Usage = Intracellular * k_usage * ( x_values[0] / Intracellular );
        RATE_OF_Intracellular_Energy = ( rate_Aerobic * 38 + rate_Anaerobic * 2 - rate_Energy_Usage ) / Intracellular;
    }


    public final double[] dy_dt_slow(double time, double[] x_values) throws Exception
    {
        this.time = time;
        this.x_values = x_values;
        final double[] dydt = new double[4];
        calculateParameters();
        calculateReactionRates();
        dydt[0] = +rate_Aerobic * 38 + rate_Anaerobic * 2 - rate_Energy_Usage; //  rate rule for of $Intracellular.Energy
        dydt[1] = -rate_Aerobic - rate_Anaerobic; //  rate rule for of $Intracellular.Glucose
        dydt[2] = +rate_Anaerobic; //  rate rule for of $Intracellular.Lactate
        dydt[3] = -rate_Aerobic * 6; //  rate rule for of $Intracellular.Oxygen
        return dydt;
    }



    @Override
    public final void init() throws Exception
    {
        CONSTRAINTS__VIOLATED = 0;
        this.simulationResultHistory.clear();
        this.simulationResultTimes.clear();
        rate_Aerobic = 0.0; // initial value of $$rate_Aerobic
        rate_Anaerobic = 0.0; // initial value of $$rate_Anaerobic
        rate_Energy_Usage = 0.0; // initial value of $$rate_Energy_Usage
        Intracellular = 1.0; // initial value of $Intracellular
        Intracellular_Lac_Secretion_Rate = 0.0; // initial value of $Intracellular.Lac_Secretion_Rate
        Intracellular_Transition_Rate = 0.0; // initial value of $Intracellular.Transition_Rate
        Intracellular_apoptosis_rate = 0.0; // initial value of $Intracellular.apoptosis_rate
        Intracellular_migration_speed = 0.0; // initial value of $Intracellular.migration_speed
        RATE_OF_Intracellular_Energy = 0.0; // initial value of RATE_OF_Intracellular_Energy
        RATE_OF_Intracellular_Glucose = 0.0; // initial value of RATE_OF_Intracellular_Glucose
        RATE_OF_Intracellular_Lactate = 0.0; // initial value of RATE_OF_Intracellular_Lactate
        RATE_OF_Intracellular_Oxygen = 0.0; // initial value of RATE_OF_Intracellular_Oxygen
        assignment = 0.0; // initial value of assignment
        energy_death_thresh = 430.0; // initial value of energy_death_thresh
        energy_prolif_thresh = 445.0; // initial value of energy_prolif_thresh
        k_aer = 0.01; // initial value of k_aer
        k_ane = 1.8E-4; // initial value of k_ane
        k_usage = 0.0023; // initial value of k_usage
        oxygen_thresh = 440.0; // initial value of oxygen_thresh
        time = 0.0; // initial value of time
        calculateInitialValues();
        initMap();
        this.isInit = true;
    }

    @Override
    protected void initMap()
    {
        variableIndex.put( "Energy", 4 );
        variableIndex.put( "Glucose", 5 );
        variableIndex.put( "Lactate", 7 );
        variableIndex.put( "Oxygen", 8 );

        variableIndex.put( "Lac_Secretion_Rate", 6 );
        variableIndex.put( "Transition_Rate", 9 );
        variableIndex.put( "apoptosis_rate", 10 );
        variableIndex.put( "migration_speed", 11 );
    }

    @Override
    public final void init(double[] initialValues, Map parameters) throws Exception
    {
        super.init( initialValues, parameters );
        this.initialValues = x_values.clone();
    }


    private final void calculateInitialValues() throws Exception
    {
        double[] x_values = this.x_values = new double[4];
        this.time = 0.0;
        x_values[0] = 450.0; //  initial value of $Intracellular.Energy
        x_values[1] = 100.0; //  initial value of $Intracellular.Glucose
        x_values[3] = 100.0; //  initial value of $Intracellular.Oxygen
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
        return new double[] {rate_Aerobic, rate_Anaerobic, rate_Energy_Usage, Intracellular, x_values[0] / Intracellular,
                x_values[1] / Intracellular, Intracellular_Lac_Secretion_Rate / Intracellular, x_values[2] / Intracellular,
                x_values[3] / Intracellular, Intracellular_Transition_Rate / Intracellular, Intracellular_apoptosis_rate / Intracellular,
                Intracellular_migration_speed / Intracellular, energy_death_thresh, energy_prolif_thresh, k_aer, k_ane, k_usage,
                oxygen_thresh, time,};
    }



    @Override
    public final void setCurrentValues(double[] values) throws Exception
    {
        CONSTRAINTS__VIOLATED = 0;
        rate_Aerobic = values[0];
        rate_Anaerobic = values[1];
        rate_Energy_Usage = values[2];
        Intracellular = values[3];
        x_values[0] = values[4] * Intracellular;
        x_values[1] = values[5] * Intracellular;
        Intracellular_Lac_Secretion_Rate = values[6] * Intracellular;
        x_values[2] = values[7] * Intracellular;
        x_values[3] = values[8] * Intracellular;
        Intracellular_Transition_Rate = values[9] * Intracellular;
        Intracellular_apoptosis_rate = values[10] * Intracellular;
        Intracellular_migration_speed = values[11] * Intracellular;
        energy_death_thresh = values[12];
        energy_prolif_thresh = values[13];
        k_aer = values[14];
        k_ane = values[15];
        k_usage = values[16];
        oxygen_thresh = values[17];
        time = values[18];
        if( time == 0 )
        {
            initialValues[0] = values[4];
            initialValues[1] = values[5];
            initialValues[2] = values[7];
            initialValues[3] = values[8];
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
        double[] z = new double[12];
        z[0] = ( NUMERIC_GT( x_values[2] / Intracellular, 0.01 ) != 0 ) ? 1 : -1; // Lac_Sec part1;
        z[1] = ( NUMERIC_GT( x_values[2] / Intracellular, 0.01 ) != 0 ) ? 1 : -1; // Lac_Sec part2;
        z[2] = ( NUMERIC_LT( x_values[0] / Intracellular, energy_death_thresh ) != 0 ) ? 1 : -1; // die part1;
        z[3] = ( NUMERIC_LT( x_values[0] / Intracellular, energy_death_thresh ) != 0 ) ? 1 : -1; // die part2;
        z[4] = ( NUMERIC_GT( x_values[0] / Intracellular, energy_prolif_thresh ) != 0 ) ? 1 : -1; // divide part1;
        z[5] = ( NUMERIC_GT( x_values[0] / Intracellular, energy_prolif_thresh ) != 0 ) ? 1 : -1; // divide part2;
        z[6] = ( NUMERIC_LT( x_values[0] / Intracellular, energy_prolif_thresh ) != 0 ) ? 1 : -1; // do_not_divide part1;
        z[7] = ( NUMERIC_LT( x_values[0] / Intracellular, energy_prolif_thresh ) != 0 ) ? 1 : -1; // do_not_divide part2;
        z[8] = ( NUMERIC_GT( x_values[0] / Intracellular, oxygen_thresh ) != 0 ) ? 1 : -1; // do_not_move part1;
        z[9] = ( NUMERIC_GT( x_values[0] / Intracellular, oxygen_thresh ) != 0 ) ? 1 : -1; // do_not_move part2;
        z[10] = ( NUMERIC_LT( x_values[0] / Intracellular, oxygen_thresh ) != 0 ) ? 1 : -1; // move part1;
        z[11] = ( NUMERIC_LT( x_values[0] / Intracellular, oxygen_thresh ) != 0 ) ? 1 : -1; // move part2;
        return z;
    }


    public final void processEvent(int i)
    {
        double[] assignments;
        double executionTime;
        switch( i )
        {
            case ( 0 ): //Lac_Sec part1
            {
                assignments = new double[2 - 1];
                executionTime = time;
                assignments[1 - 1] = 1.0E-4;
                addDelayedEvent( 0 + 1, executionTime, assignments );
                break;
            }
            case ( 1 ): //Lac_Sec part2
            {
                assignments = getNextAssignments( 1 );
                Intracellular_Lac_Secretion_Rate = assignments[0] * Intracellular;
                removeDelayedEvent( 1 );
                break;
            }
            case ( 2 ): //die part1
            {
                assignments = new double[2 - 1];
                executionTime = time;
                assignments[1 - 1] = 9.0E99;
                addDelayedEvent( 2 + 1, executionTime, assignments );
                break;
            }
            case ( 3 ): //die part2
            {
                assignments = getNextAssignments( 3 );
                Intracellular_apoptosis_rate = assignments[0] * Intracellular;
                removeDelayedEvent( 3 );
                break;
            }
            case ( 4 ): //divide part1
            {
                assignments = new double[2 - 1];
                executionTime = time;
                assignments[1 - 1] = 1.666666E-4;
                addDelayedEvent( 4 + 1, executionTime, assignments );
                break;
            }
            case ( 5 ): //divide part2
            {
                assignments = getNextAssignments( 5 );
                Intracellular_Transition_Rate = assignments[0] * Intracellular;
                removeDelayedEvent( 5 );
                break;
            }
            case ( 6 ): //do_not_divide part1
            {
                assignments = new double[2 - 1];
                executionTime = time;
                assignments[1 - 1] = 0.0;
                addDelayedEvent( 6 + 1, executionTime, assignments );
                break;
            }
            case ( 7 ): //do_not_divide part2
            {
                assignments = getNextAssignments( 7 );
                Intracellular_Transition_Rate = assignments[0] * Intracellular;
                removeDelayedEvent( 7 );
                break;
            }
            case ( 8 ): //do_not_move part1
            {
                assignments = new double[2 - 1];
                executionTime = time;
                assignments[1 - 1] = 0.0;
                addDelayedEvent( 8 + 1, executionTime, assignments );
                break;
            }
            case ( 9 ): //do_not_move part2
            {
                assignments = getNextAssignments( 9 );
                Intracellular_migration_speed = assignments[0] * Intracellular;
                removeDelayedEvent( 9 );
                break;
            }
            case ( 10 ): //move part1
            {
                assignments = new double[2 - 1];
                executionTime = time;
                assignments[1 - 1] = 10.0;
                addDelayedEvent( 10 + 1, executionTime, assignments );
                break;
            }
            case ( 11 ): //move part2
            {
                assignments = getNextAssignments( 11 );
                Intracellular_migration_speed = assignments[0] * Intracellular;
                removeDelayedEvent( 11 );
                break;
            }
            default:
                break;
        }
    }


    public final double[] getEventsPriority(double time, double[] x_values) throws Exception
    {
        calculateParameters();
        return new double[] {Double.POSITIVE_INFINITY, //Lac_Sec part1
                Double.NEGATIVE_INFINITY, //Lac_Sec part2
                Double.POSITIVE_INFINITY, //die part1
                Double.NEGATIVE_INFINITY, //die part2
                Double.POSITIVE_INFINITY, //divide part1
                Double.NEGATIVE_INFINITY, //divide part2
                Double.POSITIVE_INFINITY, //do_not_divide part1
                Double.NEGATIVE_INFINITY, //do_not_divide part2
                Double.POSITIVE_INFINITY, //do_not_move part1
                Double.NEGATIVE_INFINITY, //do_not_move part2
                Double.POSITIVE_INFINITY, //move part1
                Double.NEGATIVE_INFINITY, //move part2
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
}

