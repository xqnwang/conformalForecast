get_call_arg_with_defaults <- function(.f, call, arg) {
  defaults <- formals(.f)
  if(is.null(call[[arg]])) {
    defaults[[arg]]
  } else {
    call[[arg]]
  }
}
