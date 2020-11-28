function print_and_run_cmd(cmd::Base.AbstractCmd)

    println(cmd)

    return run(cmd)

end
