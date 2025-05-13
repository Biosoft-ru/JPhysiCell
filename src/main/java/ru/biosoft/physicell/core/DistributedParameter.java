package ru.biosoft.physicell.core;

public class DistributedParameter
{
    public String parameter;
    public Distribution distribution;


    public DistributedParameter clone()
    {
        DistributedParameter result = new DistributedParameter();
        result.distribution = distribution.clone();
        result.parameter = parameter;
        return result;
    }

    public void apply(CellDefinition cd, Model model) throws Exception
    {
        Distribution distr = distribution;
        if( distr.getName().equals( "Uniform" ) )
        {
            double min = distr.getValue( "min" );
            double max = distr.getValue( "max" );
            for( Cell cell : model.getMicroenvironment().getAgents( Cell.class ) )
            {
                if( cell.type == cd.type )
                {
                    model.signals.setSingleBehavior( cell, parameter, model.getRNG().UniformRandom( min, max ) );
                }
            }
        }
        else if( distr.getName().equals( "LogUniform" ) )
        {
            double min = distr.getValue( "min" );
            double max = distr.getValue( "max" );
            for( Cell cell : model.getMicroenvironment().getAgents( Cell.class ) )
            {
                if( cell.type == cd.type )
                {
                    model.signals.setSingleBehavior( cell, parameter, model.getRNG().LogUniformRandom( min, max ) );
                }
            }
        }
        else if( distr.getName().equals( "Normal" ) )
        {
            double min = distr.getValue( "min" );
            double max = distr.getValue( "max" );
            double mu = distr.getValue( "mu" );
            double sigma = distr.getValue( "sigma" );
            for( Cell cell : model.getMicroenvironment().getAgents( Cell.class ) )
            {
                if( cell.type == cd.type )
                {
                    model.signals.setSingleBehavior( cell, parameter, model.getRNG().NormalRandom( mu, sigma, min, max ) );
                }
            }
        }
        else if( distr.getName().equals( "LogNormal" ) )
        {
            double min = distr.getValue( "min" );
            double max = distr.getValue( "max" );
            double mu = distr.getValue( "mu" );
            double sigma = distr.getValue( "sigma" );
            for( Cell cell : model.getMicroenvironment().getAgents( Cell.class ) )
            {
                if( cell.type == cd.type )
                {
                    model.signals.setSingleBehavior( cell, parameter, model.getRNG().LogNormalRandom( mu, sigma, min, max ) );
                }
            }
        }

        else if( distr.getName().equals( "Log10Normal" ) )
        {
            double min = distr.getValue( "min" );
            double max = distr.getValue( "max" );
            double mu = distr.getValue( "mu" );
            double sigma = distr.getValue( "sigma" );
            for( Cell cell : model.getMicroenvironment().getAgents( Cell.class ) )
            {
                if( cell.type == cd.type )
                {
                    model.signals.setSingleBehavior( cell, parameter, model.getRNG().Log10NormalRandom( mu, sigma, min, max ) );
                }
            }
        }

    }
}