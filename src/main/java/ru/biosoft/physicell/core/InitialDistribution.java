package ru.biosoft.physicell.core;

public class InitialDistribution
{
    public DistributedParameter[] distributions = new DistributedParameter[0];

    public void apply(CellDefinition cDefinition, Model model) throws Exception
    {
        for( DistributedParameter sDistribution : distributions )
        {
            sDistribution.apply( cDefinition, model );
        }
    }

    public InitialDistribution clone()
    {
        InitialDistribution result = new InitialDistribution();

        result.distributions = new DistributedParameter[distributions.length];
        for( int i = 0; i < distributions.length; i++ )
            result.distributions[i] = distributions[i].clone();
        return result;
    }

}
