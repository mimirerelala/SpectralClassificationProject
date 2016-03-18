
class SetUserInRequestMiddleware:
    def process_request(self, request):
        print(request.POST.get('postname'))
        print(request.POST.get('csrfmiddlewaretoken'))
        print('Hooked in middleware')
        print(request.session.__dict__)

        request.hook = 'hooked'