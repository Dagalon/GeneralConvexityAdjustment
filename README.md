# GeneralConvexityAdjustment
This repository is associated with the paper "Convexity adjustments à la Malliavin". The link to the paper is the following:
https://arxiv.org/abs/2304.13402

We have organized the project into three parts:

ConvexityAdjustment-Lib: Here, the user can find all the necessary tools to sample the distinct processes such as SOFR, LIBOR, etc.

ConvexityAdjustmentUnitTest: We have implemented each paper example in a unit test. The file's name explains what each unit test does.

QLNet: This is Andrea Maggiulli's amazing library (https://github.com/amaggiulli/QLNet). Basically, this library is a reimplemented version of QuantLib. 
We have created a branch (https://github.com/Dagalon/QLNet.git)  to use the basic functionalities such as interest rate curves, calendars,.. all the necessary pieces to obtain the price of some interest rate derivatives.

To run the unit tests, you must follow the next steps:

1) Download the repository ConvexityAdjustment-Lib from https://github.com/Dagalon/GeneralConvexityAdjustment.git.
2) Download the branch ConvexityAdjustmentMalliavin of QLNet from https://github.com/Dagalon/QLNet.git.
3) Modify the files ConvexityAdjustment-Lib.sln and ConvexityAdjustment-Lib.csproj to specify the path of your 
QLNet's folder. This way you can add the reference to QLNet in your project.
4) Add from Nugget the library MathNet.Numeric version 5.0.0.
