function print_and_run_cmd(cmd::Cmd)

    println(cmd)

    run(cmd)

end

function print_and_run_cmd(cmd_redirect::Base.CmdRedirect)

    println(cmd_redirect)

    run(cmd_redirect)

end
