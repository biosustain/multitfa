from functools import wraps 
""" Adopted from https://medium.com/@mgarod/dynamically-add-a-method-to-a-class-in-python-c49204b85bd6

Defining a decorator that accepts the method for a class

"""
def add_method(cls):
    def decorator(func):
        @wraps(func) 
        def wrapper(self, *args, **kwargs): 
            return func(*args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        # Note we are not binding func, but wrapper which accepts self but does exactly the same as func
        return func # returning func means func can still be used normally
    return decorator