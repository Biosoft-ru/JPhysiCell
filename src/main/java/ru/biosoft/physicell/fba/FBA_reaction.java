package ru.biosoft.physicell.fba;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class FBA_reaction
{
    String id;
    String name;

    double lowerBound;
    double upperBound;
    double objectiveCoefficient;
    double fluxValue;

    Map<FBA_metabolite, Double> metabolites;
    Map<String, FBA_metabolite> idMetaboliteMap;

    FBA_reaction(String id)
    {
        this.id = id;
        this.name = "";
        this.lowerBound = 0;
        this.upperBound = 1000;
        this.objectiveCoefficient = 0;
        this.fluxValue = 0;
    }


    String getId()
    {
        return this.id;
    }

    void setName(String name)
    {
        this.name = name;
    }

    String getName()
    {
        return this.name;
    }

    void setLowerBound(double lowerBound)
    {
        this.lowerBound = lowerBound;
    }

    double getLowerBound()
    {
        return this.lowerBound;
    }

    void setUpperBound(double upperBound)
    {
        this.upperBound = upperBound;
    }

    double getUpperBound()
    {
        return this.upperBound;
    }

    void setObjectiveCoefficient(double ojectiveCoefficient)
    {
        this.objectiveCoefficient = ojectiveCoefficient;
    }

    double getObjectiveCoefficient()
    {
        return this.objectiveCoefficient;
    }

    void setFluxValue(double fluxValue)
    {
        this.fluxValue = fluxValue;
    }

    public double getFluxValue()
    {
        return this.fluxValue;
    }

    int getNumberOfMetabolites()
    {
        return this.metabolites.size();
    }

    Map<FBA_metabolite, Double> getMetabolites()
    {
        return this.metabolites;
    }

    boolean reversible()
    {
        return lowerBound < 0;
    }

    boolean hasMetabolite(String mId)
    {
        return idMetaboliteMap.containsKey( mId );
    }

    void addMetabolite(FBA_metabolite met, double stoich)
    {
        if( this.hasMetabolite( met.getId() ) )
        {
            double val = metabolites.get( met );
            this.metabolites.put( met, val + stoich );
        }
        else
        {
            this.idMetaboliteMap.put( met.getId(), met );
            this.metabolites.put( met, stoich );
        }
    }

    List<String> getReactants()
    {
        List<String> reactants = new ArrayList<>();
        for( Entry<FBA_metabolite, Double> entry : metabolites.entrySet() )
        {
            FBA_metabolite met = entry.getKey();
            double sotich = entry.getValue();
            if( sotich < 0 )
            {
                String mId = met.getId();
                reactants.add( mId );
            }
        }

        return reactants;
    }

    List<String> getProducts()
    {
        List<String> products = new ArrayList<>();
        //        for( auto itr = this.metabolites.begin(); itr != this.metabolites.end(); ++itr )
        //        {
        for( Entry<FBA_metabolite, Double> entry : metabolites.entrySet() )
        {
            FBA_metabolite met = entry.getKey();//.first;
            double sotich = entry.getValue();//.second;
            if( sotich > 0 )
            {
                String mId = met.getId();
                products.add( mId );
            }
        }
        return products;
    }

    double getStoichCoefficient(String mId)
    {
        double coefficient = 0;
        if( hasMetabolite( mId ) )
        {
            FBA_metabolite met = this.idMetaboliteMap.get( mId );
            coefficient = this.metabolites.get( met );
        }
        return coefficient;
    }

    String getReactionString()
    {
        List<String> compounds;
        String reactionString = "";

        compounds = this.getReactants();
        for( int i = 0; i < compounds.size(); i++ )
        {
            String mId = compounds.get( i );
            FBA_metabolite met = this.idMetaboliteMap.get( mId );
            // Change the sign (-) of the coeeficient for printing
            double coeff = -1. * this.getStoichCoefficient( mId );
            if( (int)coeff == coeff )
                reactionString += (int)coeff;
            else
                reactionString += coeff;

            //reactionString += " " + met.getName();
            reactionString += " " + met.getId();

            if( i < compounds.size() - 1 )
                reactionString += " + ";
        }
        if( reversible() )
            reactionString += " <==> ";
        else
            reactionString += " -. ";

        compounds = getProducts();
        for( int i = 0; i < compounds.size(); i++ )
        {
            String mId = compounds.get( i );
            FBA_metabolite met = this.idMetaboliteMap.get( mId );
            double coeff = getStoichCoefficient( mId );
            if( (int)coeff == coeff )
                reactionString += (int)coeff;
            else
                reactionString += coeff;

            reactionString += " " + met.getId();

            if( i < compounds.size() - 1 )
                reactionString += " + ";
        }
        return reactionString;
    }
}