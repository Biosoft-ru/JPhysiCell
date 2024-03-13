package ru.biosoft.physicell.core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class HypothesisRuleset
{
    Map<String, HypothesisRule> rulesMap = new HashMap<>();
    String type;
    CellDefinition cd;
    List<HypothesisRule> rules = new ArrayList<>();

    HypothesisRuleset()
    {
        type = "none";
        cd = null;
        rules = new ArrayList<>();
        rulesMap.clear();
    }

    public String display()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "\nBehavioral rules for cell type " + type + ":" );
        sb.append( "\n===================================================" );
        for( HypothesisRule rule : rules )
            sb.append( rule.reducedDisplay() );
        sb.append( "\n" );
        return sb.toString();
    }

    String detailedDisplay()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "\nBehavioral rules for cell type " + type + ":" );
        sb.append( "\n===================================================" );
        for( HypothesisRule rule : rules )
            sb.append( rule.detailedDisplay() );
        sb.append( "\n" );
        return sb.toString();
    }

    void sync(CellDefinition cd)
    {
        this.cd = cd;
        type = cd.name;
        for( HypothesisRule rule : rules )
            rule.sync( cd );
    }

    HypothesisRule addBehavior(String behavior, double minBehavior, double maxBehavior) throws Exception
    {
        // check: is this a valid signal? (is it in the dictionary?)
        if( SignalBehavior.findBehaviorIndex( behavior ) < 0 )
            throw new Exception( "Warning! Attempted to add behavior " + behavior + " which is not in the dictionary. "
                    + "Either fix your model or add the missing behavior to the simulation." );

        // first, check. Is there already a ruleset? 
        HypothesisRule rule = rulesMap.get( behavior );

        // if not, add it 
        if( rule == null )
        {
            HypothesisRule hr = new HypothesisRule();
            hr.behavior = behavior;
            hr.sync( cd );
            hr.minValue = minBehavior;
            hr.maxValue = maxBehavior;
            rules.add( hr );
            rulesMap.put( behavior, hr );
            return hr;
        }

        // otherwise, edit it 
        //        HypothesisRule pHR = search.second; 

        /*
            // March 28 2023 fix  : let's not overwrite eixsting values
        pHR.min_value = min_behavior; 
        pHR.max_value = max_behavior; 
        */

        return rule;
    }

    HypothesisRule addBehavior(String behavior) throws Exception
    {
        double minBehavior = 0.1;
        double maxBehavior = 1.0;
        return addBehavior( behavior, minBehavior, maxBehavior );
    }

    void sync(String name)
    {
        sync( CellDefinition.getCellDefinition( name ) );
    }

    HypothesisRule findBehavior(String name)
    {
        return rulesMap.get( name );
    }

    HypothesisRule get(String name)
    {
        return findBehavior( name );
    }

    void apply(Cell cell) throws Exception
    {
        for( HypothesisRule rule : rules )
            rule.apply( cell );
    }
}
