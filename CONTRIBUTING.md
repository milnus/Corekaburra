# Welcome to the Corekaburra contributing guide

If you are looking to contribute to Corekaburra, thank you! Any contribution you make is greatly appreciated and will be acknowledged. If you are looking to report an issuse, implement a feature, or in general discuss the workings of Corekaburra you can find ways of doing this below. 

Read our [Code of Conduct](./CODE_OF_CONDUCT.md) to keep our community approachable and respectable.

## Raise an issue
If you have observed a problem with Corekaburra there are a few things to check:
1. Do you have the latest version installed? - If not make sure to test the latest version, as things may have been fixed.
2. Search if a similar issue has been raised on the [issues tab](https://github.com/milnus/Corekaburra/issues) and if there are any solutions (remember to also look at closed issues). - If a similar issue has not had a solution, post a comment in the issue that you have experienced a similar problem.
3. Try to pin down the problem to the best of your ability. The more precise you can be the easier and faster we can get to the root cause and find a fix. If you are able to construct a small example that can reproduce your issue we greatly appreciate it!

When you raise an issue on the [issue tab of Corekaburra](https://github.com/milnus/Corekaburra/issues) we will try to get around to you as quickly as possible, but we may have other pressing tasks to tend to. 

## Contribute code
If you have an issue and/or a feature you miss in Corekaburra, and feel like you can/have time to code it yourself, here are some guidelines:

1. Keep to the dependencies. We do not want to add anymore dependencies to Corekaburra, therefore do not write code that makes use of packages outside of standard python packages, Networkx, biopython, gffutils, and numpy.
2. Describe your code. Use doc strings in the beginning of each function marked by triple quotes ("""doc-string""") described the intention of the function, input arguments and returned objects.
3. No tests required. We do not expect you to write any unit- or functional tests for you code. We would like if you break up your code into functions that makes it easier to produce unittest for.
5. Make a draft pull-request. By making a draft pull-request to the main branch you can get the [Github action](https://github.com/milnus/Corekaburra/actions) for Corekaburra to run the unit and functional tests, when you commit your changes. It is not essential to pass the tests, as we can fix things up for you before merging in your code. 
4. Make a pull request. When you have implemented your code, make a pull request with a short description of the aim of your contribution and what major things that have been changed. We will review the changes and merge when approved.

## Discuss workings
If you have thoughts to how Corekaburra handles certain things raise an issue on the [issue tab](https://github.com/milnus/Corekaburra/issues). Be mindful, respectful and inviting of discussion. We learn all the time and would like to improve the tool, so we are keen to learn and discuss subjects, but only in a respectful tone.
