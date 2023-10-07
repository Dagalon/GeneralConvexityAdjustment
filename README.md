# GeneralConvexityAdjustment
This repository is associated to the paper "Convexity adjustments à la Malliavin". The link of the paper is the next:

https://arxiv.org/abs/2304.13402

The project has been organized in three parts:

ConvexityAdjustment-Lib: Where the user can find all tools that we will need to sampling the distinct processes SOFR, LIBOR,...

ConvexityAdjustmentUnitTest: Each the example of the paper has been implemented in a unit test. We think, that The name 
of the file explains what each unit test does.

QLNet: This is the amazing library of Andrea Maggiulli (https://github.com/amaggiulli/QLNet). Basically, this library is reimplemented version of QuantLib. 
We have created a branch (https://github.com/Dagalon/QLNet.git) of QLNet in order to be able to use the basic functionalities as interest rate curve, calendars,..
that you can need when you try to get the price of some interest rate derivative.

In order to run the unit tests, you must follow the next steps:

1) Download the repository ConvexityAdjustment-Lib from https://github.com/Dagalon/GeneralConvexityAdjustment.git.
2) Download the branch ConvexityAdjustmentMalliavin of QLNet from https://github.com/Dagalon/QLNet.git.
3) Modify the files ConvexityAdjustment-Lib.sln and ConvexityAdjustment-Lib.csproj in order to specify the path of your 
QLNet's folder, of this way you add the reference to QLNet in your project.
4) Add from Nugget the library MathNet.Numeric version 5.0.0.

