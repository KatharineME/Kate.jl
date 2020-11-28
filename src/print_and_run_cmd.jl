function print_and_run_cmd(cmd::Base.AbstractCmd)

    println(cmd)

    run(cmd)

end
