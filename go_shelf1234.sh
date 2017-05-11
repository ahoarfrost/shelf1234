#go_shelf1234.sh

mkdir figures

cd FLA
RScript FlaRates_shelf1234.R
RScript FlaRatesSpecialCases_shelf1234.R
RScript FlaMaxesFactors_shelf1234.R
RSript figuresFLA_shelf1234.R

cd ../plate
mkdir flvstime-plots
RScript rawToRates_shelf1234.R
RScript furtherProcessing_shelf1234.R
RScript figuresPlate_shelf1234.R
