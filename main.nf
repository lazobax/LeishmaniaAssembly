include {Test} from './modules/test'
include {Test2} from './modules/test2'

workflow{

    Test("a")
    Test2(Test.out, "b")


}